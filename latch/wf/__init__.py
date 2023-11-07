"""Latch wrapper of ArchR plotEmbedding function.
"""

import subprocess
import glob

from pathlib import Path
from dataclasses import dataclass
from typing import List

from latch import large_task, small_task, workflow
from latch.registry.table import Table
from latch.types import (
    LatchAuthor,
    LatchDir,
    LatchMetadata,
    LatchParameter,
)


@dataclass
class Run:
    run_name: str
    archrObj: LatchDir


@dataclass
class Project:
    project_id: str
    archrObj: LatchDir


@dataclass
class ShinyProject:
    shiny_dir: LatchDir
    run_names: List[str]


def initialize_runs(
    projects: List[Project],
    project_table_id: str
) -> List[Run]:

    runs = []
    project_table = Table(project_table_id)
    try:
        for p in projects:
            for page in project_table.list_records():
                for p_id, record in page.items():
                    p_id = record.get_name()
                    if p_id == p.project_id:
                        project = record.get_values()
                        try:
                            if len(project['Runs']) > 0:
                                for project_run in project['Runs']:
                                    run_info = project_run.get_values()
                                    run_id = project_run.get_name()
                                    try:
                                        run_archr = run_info[
                                            'archrproject_outs'
                                        ]
                                        runs.append(Run(run_id, run_archr))
                                    except:
                                        print(
                                            f"Data missing for run: {run_id}"
                                        )
                        except:
                            break
        return runs
    except Exception as err:
        print(f"Unexpected {err=}, {type(err)=}")
        return


@large_task
def runScript(
    projects: List[Project],
    project_table_id: str,
    output_dir: LatchDir,
    project: str = "test",
    groupBy: str = "Clusters",
) -> ShinyProject:

    runs = initialize_runs(projects, project_table_id)

    local_archr_dir = projects[0].archrObj.local_path

    r_names_list = []
    for r in runs:
        r_names_list.append(r.run_name)

    run_ids = [
        Path(i).stem for i in glob.glob(
            f"{local_archr_dir}/**/*.arrow", recursive=True
        )
    ]

    subprocess.run(
        ["Rscript", "/root/wf/runShiny.R", local_archr_dir, project, groupBy]
    )

    local_output_dir = str(Path("/root/").resolve())

    remote_path = output_dir.remote_path
    if remote_path[-1] != "/":
        remote_path += "/"

    return ShinyProject(
        shiny_dir=LatchDir(local_output_dir, remote_path),
        run_names=r_names_list
    )


@small_task
def upload_to_registry(
    projects: List[Project],
    shiny_project: ShinyProject,
    run_table_id: str = "761",
    project_table_id: str = "917",
):

    runs = initialize_runs(projects, project_table_id)
    run_table = Table(run_table_id)
    project_table = Table(project_table_id)

    try:
        with run_table.update() as updater:
            for run in shiny_project.run_names:
                updater.upsert_record(
                    run,
                    atlasshiny_outs=shiny_project.shiny_dir
                )

        with project_table.update() as updater:
            for p in projects:
                updater.upsert_record(
                    p.project_id,
                    atlasshiny_outs=shiny_project.shiny_dir
                )

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

        "projects": LatchParameter(
            display_name="ArchR Objects",
            description="Select projects for which the ArchRObject will be \
                        processed.",
            samplesheet=True,
            batch_table_column=True,
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
            display_name="Output Directory",
            description="Where to save the plots?."
        ),
        "run_table_id": LatchParameter(
            display_name="Registry Table ID",
            description="Table in Registry to save the Shiny output to",
        ),
        'project_table_id': LatchParameter(
            display_name='The ID of the SOWs Registry table',
            description='The Shiny project will be inserted into the SOW \
                        table for the corresponding runs.'
        )
    },
)


@workflow(metadata)
def shinyArchr_wf(
    projects: List[Project],
    output_dir: LatchDir,
    project: str = "test",
    groupBy: str = "Clusters",
    run_table_id: str = "761",
    project_table_id: str = "917"
) -> ShinyProject:
    """is a full-featured software suite for the analysis of single-cell
    chromatin accessibility data.

    atlasShiny
    ----

    `atlasShiny` is a full-featured application for the exploring of ArchR
    output data.
    """

    shiny_project = runScript(
        output_dir=output_dir,
        project=project,
        groupBy=groupBy,
        projects=projects
    )

    upload_to_registry(
        projects=projects,
        shiny_project=shiny_project,
        run_table_id=run_table_id,
        project_table_id=project_table_id
    )

    return shiny_project


# if __name__ == "__main__":

#     shinyArchr_wf(
#         projects=[Project("sample_project\n\n", LatchDir())],
#         output_dir=LatchDir("latch://13502.account/"),
#         project="D1266_w_chromap_frag",
#         groupBy="Clusters",
#         run_table_id="761",
#         project_table_id="917"
#     )
