"""
Latch wrapper of ArchR plotEmbedding function.
"""


import subprocess
import glob
from pathlib import Path

from flytekit import LaunchPlan, task, workflow
from latch.types import LatchDir
from latch.types import LatchFile
from latch.registry.table import Table
from dataclasses import dataclass
from dataclasses_json import dataclass_json
from enum import Enum
from typing import List
from latch import large_task, small_task, workflow

from latch.types import (
    LatchAuthor,
    LatchDir,
    LatchFile,
    LatchMetadata,
    LatchParameter,
    LatchRule,
)

@dataclass
class RunName:
    run_name: str
    archrObj: LatchDir

# @dataclass
# class ArchrObject:
#     archrObj: LatchDir


@dataclass
class ShinyProject:
    run_ids: List[str]
    shiny_dir: LatchDir
    run_names: List[str]


@large_task
def runScript(
    runNames: List[RunName],
    output_dir: LatchDir,
    project: str = "test",
    groupBy: str = "Clusters",
) -> ShinyProject:
    local_archr_dir = runNames[0].archrObj.local_path

    print(local_archr_dir)
    print("RUN NAMES")
    r_names_list = []
    for r in runNames:
        print(r.run_name)
        r_names_list.append(r.run_name)
    print(runNames)

    run_ids = [
        Path(i).stem for i in glob.glob(f"{local_archr_dir}/**/*.arrow", recursive=True)
    ]

    subprocess.run(
        ["Rscript", "/root/wf/runShiny.R", local_archr_dir, project, groupBy]
    )

    local_output_dir = str(Path(f"/root/").resolve())

    remote_path = output_dir.remote_path
    if remote_path[-1] != "/":
        remote_path += "/"

    return ShinyProject(
        run_ids=run_ids, shiny_dir=LatchDir(local_output_dir, remote_path),run_names=r_names_list
    )

@small_task
def upload_to_registry(
    shiny_project: ShinyProject,
    run_table_id: str = "761",
    project_table_id: str = "779",
):
    run_table = Table(run_table_id)
    project_table = Table(project_table_id)

    try:
        with run_table.update() as updater:
            for run in shiny_project.run_names:
                updater.upsert_record(run, atlasshiny_outs=shiny_project.shiny_dir)

        with project_table.update() as updater:
            for run in shiny_project.run_names:
                updater.upsert_record(run, atlasshiny_outs=shiny_project.shiny_dir)

    except Exception as err:
        print(f"Unexpected {err=}, {type(err)=}")
        return
    finally:
        return

metadata = LatchMetadata(
    display_name="atlasShiny",
    author=LatchAuthor(
        name="Noori",
        email="noorisotude@gmail.com",
        github="github.com/atlasxomics",
    ),
    repository="https://github.com/atlasxomics/ShinyCell",
    license="MIT",
    parameters={

        "runNames": LatchParameter(
            display_name="ArchR Objects",
            description="Select runs",
            samplesheet=True,
        ),
        "project": LatchParameter(
            display_name="Project Name",
            description="Specify a name for the output folder.",
        ),
        "groupBy": LatchParameter(
            display_name="Group By",
            description="A string that indicates how cells should be grouped.",
        ),
        "output_dir": LatchParameter(
            display_name="Output Directory", description="Where to save the plots?."
        ),
        "run_table_id": LatchParameter(
            display_name="Registry Table ID",
            description="Table in Registry to save the Shiny output to",
        ),
        'project_table_id': LatchParameter(
            display_name='The ID of the SOWs Registry table',
            description='The Shiny project will be inserted into the SOW table for the corresponding runs.'
        )
    },
)

@workflow(metadata)
def shinyArchr_wf(
    runNames: List[RunName],
    output_dir: LatchDir,
    project: str = "test",
    groupBy: str = "Clusters",
    run_table_id: str = "761",
    project_table_id: str = "779"
) -> ShinyProject:
    """is a full-featured software suite for the analysis of single-cell chromatin accessibility data.

    atlasShiny
    ----

    `atlasShiny` is a full-featured application for the exploring of ArchR output data.
    """

    shiny_project = runScript(
        output_dir=output_dir, project=project, groupBy=groupBy, runNames=runNames
    )

    upload_to_registry(shiny_project=shiny_project, run_table_id=run_table_id, project_table_id=project_table_id)

    return shiny_project

if __name__ == "__main__":

    shinyArchr_wf(
        output_dir=LatchDir("latch://13502.account/rshinyA_outs"),
        project="D1266_w_chromap_frag",
        groupBy="Clusters",
        run_table_id="761",
        project_table_id="779"
    )