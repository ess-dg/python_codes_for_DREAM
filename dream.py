import pandas as pd
import numpy as np
from nexusutils.nexusbuilder import NexusBuilder

"""
Generates mesh geometry for DREAM Endcap detector from information from a GEANT4 simulation
"""


def find_voxel_vertices(
        dz: float,
        theta: float,
        phi: float,
        dy1: float,
        dx1: float,
        dx2: float,
        alp1: float,
        dy2: float,
        dx3: float,
        dx4: float,
        alp2: float,
) -> (float, float, float, float, float, float, float, float):
    """
    Ported from GEANT4
    http://www.apc.univ-paris7.fr/~franco/g4doxy/html/G4Trap_8cc-source.html
    http://geant4-userdoc.web.cern.ch/geant4-userdoc/UsersGuides/ForApplicationDeveloper/html/Detector/Geometry/geomSolids.html
    """
    ttheta_cphi = np.tan(theta) * np.cos(phi)
    ttheta_sphi = np.tan(theta) * np.sin(phi)
    talpha1 = np.tan(alp1)
    talpha2 = np.tan(alp2)

    pt_0 = np.array(
        (
            -dz * ttheta_cphi - dy1 * talpha1 - dx1,
            -dz * ttheta_sphi - dy1,
            -dz,
        )
    )
    pt_1 = np.array(
        (
            -dz * ttheta_cphi - dy1 * talpha1 + dx1,
            -dz * ttheta_sphi - dy1,
            -dz,
        )
    )
    pt_2 = np.array(
        (
            -dz * ttheta_cphi + dy1 * talpha1 - dx2,
            -dz * ttheta_sphi + dy1,
            -dz,
        )
    )
    pt_3 = np.array(
        (
            -dz * ttheta_cphi + dy1 * talpha1 + dx2,
            -dz * ttheta_sphi + dy1,
            -dz,
        )
    )
    pt_4 = np.array(
        (
            +dz * ttheta_cphi - dy2 * talpha2 - dx3,
            +dz * ttheta_sphi - dy2,
            +dz,
        )
    )
    pt_5 = np.array(
        (
            +dz * ttheta_cphi - dy2 * talpha2 + dx3,
            +dz * ttheta_sphi - dy2,
            +dz,
        )
    )
    pt_6 = np.array(
        (
            +dz * ttheta_cphi + dy2 * talpha2 - dx4,
            +dz * ttheta_sphi + dy2,
            +dz,
        )
    )
    pt_7 = np.array(
        (
            +dz * ttheta_cphi + dy2 * talpha2 + dx4,
            +dz * ttheta_sphi + dy2,
            +dz,
        )
    )
    return pt_0, pt_1, pt_2, pt_3, pt_4, pt_5, pt_6, pt_7


def create_winding_order(
        number_of_voxels: int,
        vertices_in_voxel: int,
        vertices_in_each_face: int,
        vertex_start_index: int,
) -> np.ndarray:
    index_0 = []
    index_1 = []
    index_2 = []
    index_3 = []
    for voxel in range(number_of_voxels):
        start_index = (voxel * vertices_in_voxel) + vertex_start_index
        index_0.extend(
            [
                start_index,
                start_index,
                start_index,
                start_index + 1,
                start_index + 2,
                start_index + 4,
            ]
        )
        index_1.extend(
            [
                start_index + 2,
                start_index + 4,
                start_index + 1,
                start_index + 3,
                start_index + 6,
                start_index + 5,
            ]
        )
        index_2.extend(
            [
                start_index + 3,
                start_index + 6,
                start_index + 5,
                start_index + 7,
                start_index + 7,
                start_index + 7,
            ]
        )
        index_3.extend(
            [
                start_index + 1,
                start_index + 2,
                start_index + 4,
                start_index + 5,
                start_index + 3,
                start_index + 6,
            ]
        )

    data = np.column_stack(
        (
            vertices_in_each_face,
            index_0,
            index_1,
            index_2,
            index_3,
        )
    ).astype(np.int32)
    return data


def write_to_off_file(
        filename: str,
        number_of_vertices: int,
        number_of_faces: int,
        vertices: np.ndarray,
        voxels: np.ndarray,
):
    """
    Write mesh geometry to a file in the OFF format
    https://en.wikipedia.org/wiki/OFF_(file_format)
    """
    with open(filename, "w") as f:
        f.writelines(
            (
                "OFF\n",
                "# DREAM Detectors\n",
                f"{number_of_vertices} {number_of_faces} 0\n",
            )
        )
    with open(filename, "a") as f:
        pd.DataFrame(vertices).to_csv(f, sep=" ", header=None, index=False)
    with open(filename, "a") as f:
        pd.DataFrame(voxels).to_csv(f, sep=" ", header=None, index=False)


def rotate_around_x(angle_degrees: float, vertex: np.ndarray) -> np.ndarray:
    angle = np.deg2rad(angle_degrees)
    rotation_matrix = np.array(
        [
            [1, 0, 0],
            [0, np.cos(angle), -np.sin(angle)],
            [0, np.sin(angle), np.cos(angle)],
        ]
    )
    return rotation_matrix.dot(vertex)


def rotate_around_y(angle_degrees: float, vertex: np.ndarray) -> np.ndarray:
    angle = np.deg2rad(angle_degrees)
    rotation_matrix = np.array(
        [
            [np.cos(angle), 0, np.sin(angle)],
            [0, 1, 0],
            [-np.sin(angle), 0, np.cos(angle)],
        ]
    )
    return rotation_matrix.dot(vertex)


def rotate_around_z(angle_degrees: float, vertex: np.ndarray) -> np.ndarray:
    angle = np.deg2rad(angle_degrees)
    rotation_matrix = np.array(
        [
            [np.cos(angle), -np.sin(angle), 0],
            [np.sin(angle), np.cos(angle), 0],
            [0, 0, 1],
        ]
    )
    return rotation_matrix.dot(vertex)

# Irina: added system_id to the list, identifier for the DREAM detector type
# system_id                  Dream detector
#    3                    SUMO3, EndCap Backward
#    4                    SUMO4, EndCap Backward
#    5                    SUMO5, EndCap Backward
#    6                    SUMO6, EndCap Backward
#    7                    Mantle
#    8                    High-Resolution
#    9                    SANS
#    13                   SUMO3, EndCap Forward
#    14                   SUMO4, EndCap Forward
#    15                   SUMO5, EndCap Forward
#    16                   SUMO6, EndCap Forward

def write_to_nexus_file(
        filename: str,
        system_id: np.ndarray,
        number_of_vertices: int,
        vertices: np.ndarray,
        voxels: np.ndarray,
        detector_ids: np.ndarray,
        wire_id: np.ndarray,
        strip_id: np.ndarray,
        x_offsets: np.ndarray,
        y_offsets: np.ndarray,
        z_offsets: np.ndarray,
):
    vertices_in_face = 4
    faces = np.arange(0, number_of_vertices, vertices_in_face)

    with NexusBuilder(
            filename, compress_type="gzip", compress_opts=1, nx_entry_name="entry"
    ) as builder:
        instrument_group = builder.add_nx_group(builder.root, "DREAM", "NXinstrument")
        builder.add_dataset(instrument_group, "name", "DREAM")

        detector_group = builder.add_nx_group(
            instrument_group, "all_detectors", "NXdetector"
        )

        # system_id_group = builder.add_nx_group(detector_group, 'system_ID', "NXlog")
        # builder.add_dataset(system_id_group, "sumo_number", system_id.astype(np.int32))

        builder.add_dataset(detector_group, "sumo_number", system_id.astype(np.int32))

        builder.add_dataset(detector_group, "detector_number", np.unique(detector_ids[:, 1]).astype(np.int32))
        builder.add_dataset(detector_group, "wire_number", wire_id.astype(np.int32))
        builder.add_dataset(detector_group, "strip_number", strip_id.astype(np.int32))

        # Record voxel centre positions
        builder.add_dataset(detector_group, "x_pixel_offset", x_offsets, {"units": "mm"})
        builder.add_dataset(detector_group, "y_pixel_offset", y_offsets, {"units": "mm"})
        builder.add_dataset(detector_group, "z_pixel_offset", z_offsets, {"units": "mm"})

        builder.add_fake_event_data(1, 1000)

        shape_group = builder.add_nx_group(
            detector_group, "detector_shape", "NXoff_geometry"
        )

        builder.add_dataset(shape_group, "vertices", vertices.astype(np.float64))
        builder.add_dataset(
            shape_group, "winding_order", voxels.flatten().astype(np.int32)
        )
        builder.add_dataset(shape_group, "faces", faces.astype(np.int32))
        builder.add_dataset(
            shape_group, "detector_faces", detector_ids.astype(np.int32)
        )


def create_sector(
        geant_df: pd.DataFrame,
        max_vertex_index: int,
        max_face_index: int,
) -> (np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray):
    number_of_voxels = len(df.index)
    vertices_in_voxel = 8
    faces_in_voxel = 6
    number_of_vertices = vertices_in_voxel * number_of_voxels
    number_of_faces = faces_in_voxel * number_of_voxels

    detector_number = np.zeros(number_of_voxels)

    wire_number = np.zeros(number_of_voxels)
    strip_number = np.zeros(number_of_voxels)

    x_coords = np.zeros(number_of_vertices)
    y_coords = np.zeros(number_of_vertices)
    z_coords = np.zeros(number_of_vertices)

    x_centre_coords = np.zeros(number_of_voxels)
    y_centre_coords = np.zeros(number_of_voxels)
    z_centre_coords = np.zeros(number_of_voxels)

    voxel_ids = np.zeros((number_of_faces, 2))

    max_voxel_index = max_face_index / faces_in_voxel

    for voxel in range(number_of_voxels):
        # Map each face in the voxel to the voxel ID
        for face_number_in_voxel in range(faces_in_voxel):
            face = voxel * faces_in_voxel + face_number_in_voxel + max_face_index
            voxel_ids[voxel * faces_in_voxel + face_number_in_voxel, 0] = face
            voxel_ids[voxel * faces_in_voxel + face_number_in_voxel, 1] = (
                    voxel + max_voxel_index
            )

        voxel_vertices = find_voxel_vertices(
            geant_df["z"][voxel] / 2,
            geant_df["dtheta"][voxel],
            0.0,
            geant_df["y1"][voxel] / 2,
            geant_df["x1"][voxel] / 2,
            geant_df["x1"][voxel] / 2,
            0.0,
            geant_df["y2"][voxel] / 2,
            geant_df["x2"][voxel] / 2,
            geant_df["x2"][voxel] / 2,
            0.0,
        )

        voxel_position = np.array(
            [
                geant_df["x_centre"][voxel],
                geant_df["y_centre"][voxel],
                geant_df["z_centre"][voxel],
            ]
        )

        voxel_rotation = np.array(
            [
                [geant_df["rxx"][voxel], geant_df["rxy"][voxel], geant_df["rxz"][voxel]],
                [geant_df["ryx"][voxel], geant_df["ryy"][voxel], geant_df["ryz"][voxel]],
                [geant_df["rzx"][voxel], geant_df["rzy"][voxel], geant_df["rzz"][voxel]],
            ]
        )

        detector_number[voxel] = geant_df["sumo"][voxel]

        wire_number[voxel] = geant_df["wire"][voxel]
        strip_number[voxel] = geant_df["strip"][voxel]

        for vert_number, vertex in enumerate(voxel_vertices):
            # Rotate voxel to final position, added by IS
            # detector_number = 7 is Mantle
            if detector_number[voxel] == 7:
                vertex = rotate_around_y(voxel_rotation[0, 0], vertex)
                vertex = rotate_around_z(voxel_rotation[0, 1], vertex)
                vertex = rotate_around_z(voxel_rotation[0, 2], vertex)
            # detector_number = 8 is High-Resolution detector amd 9 for SANS detector
            elif detector_number[voxel] in [8, 9]:
                vertex = rotate_around_z(voxel_rotation[0, 0], vertex)
                vertex = rotate_around_x(voxel_rotation[0, 1], vertex)
                vertex = rotate_around_y(voxel_rotation[0, 2], vertex)
                vertex = rotate_around_z(voxel_rotation[1, 0], vertex)
            # EndCap Backward (det_num = 3-6) and Forward (13-16)
            else:
                # vertex = rotate_around_y(0, vertex)
                vertex = rotate_around_y(voxel_rotation[0, 0], vertex)
                vertex = rotate_around_z(voxel_rotation[0, 1], vertex)
                vertex = rotate_around_y(voxel_rotation[0, 2], vertex)
                vertex = rotate_around_x(voxel_rotation[1, 0], vertex)
                vertex = rotate_around_z(voxel_rotation[1, 1], vertex)
                vertex = rotate_around_y(voxel_rotation[1, 2], vertex)
                vertex = rotate_around_z(voxel_rotation[2, 0], vertex)

            vertex += voxel_position

            x_coords[voxel * vertices_in_voxel + vert_number] = vertex[0]
            y_coords[voxel * vertices_in_voxel + vert_number] = vertex[1]
            z_coords[voxel * vertices_in_voxel + vert_number] = vertex[2]

        x_centre_coords[voxel] = np.mean(
            x_coords[
              voxel * vertices_in_voxel: voxel * vertices_in_voxel
                                         + vertices_in_voxel
            ]
        )
        y_centre_coords[voxel] = np.mean(
            y_coords[
              voxel * vertices_in_voxel: voxel * vertices_in_voxel
                                         + vertices_in_voxel
            ]
        )
        z_centre_coords[voxel] = np.mean(
            z_coords[
              voxel * vertices_in_voxel: voxel * vertices_in_voxel
                                         + vertices_in_voxel
            ]
        )

    vertex_coords = np.column_stack((x_coords, y_coords, z_coords))

    # Vertices making up each face of each voxel

    vertices_in_each_face = 4 * np.ones(number_of_faces)

    faces = create_winding_order(
        number_of_voxels, vertices_in_voxel, vertices_in_each_face, max_vertex_index
    )
    return (
        detector_number,
        wire_number,
        strip_number,
        vertex_coords,
        faces,
        voxel_ids,
        x_centre_coords,
        y_centre_coords,
        z_centre_coords,
    )


if __name__ == "__main__":
    df = pd.read_csv(
        "DREAMAll_voxels.txt",
        delim_whitespace=True,
        header=None
    )
    df.columns = [
        "sumo",
        "sect",
        "module",
        "seg",
        "wire",
        "strip",
        "counter",
        "x_centre",
        "y_centre",
        "z_centre",
        "dtheta",
        "x1",
        "x2",
        "y1",
        "y2",
        "z",
        "rxx",
        "rxy",
        "rxz",
        "ryx",
        "ryy",
        "ryz",
        "rzx",
        "rzy",
        "rzz",
    ]

    faces_in_voxel = 6

    det_number = None
    w_number = None
    st_number = None
    total_vertices = None
    total_faces = None
    total_ids = None
    x_offsets_total = None
    y_offsets_total = None
    z_offsets_total = None
    max_vertex_index = 0
    max_face_index = 0

    # TODO start and stop angle are inferred from diagrams, need to check
    sector_number, wire_number, strip_number, sector_vertices, sector_faces, sector_ids, x_offsets, y_offsets, z_offsets = create_sector(
        df,
        max_vertex_index,
        max_face_index,
    )

    if total_vertices is None:
        det_number = sector_number
        w_number = wire_number
        st_number = strip_number
        total_vertices = sector_vertices
        total_faces = sector_faces
        total_ids = sector_ids
        x_offsets_total = x_offsets
        y_offsets_total = y_offsets
        z_offsets_total = z_offsets
    else:
        det_number = np.vstack((det_number, sector_number))
        w_number = np.vstack((w_number, wire_number))
        st_number = np.vstack((st_number, strip_number))
        total_vertices = np.vstack((total_vertices, sector_vertices))
        total_faces = np.vstack((total_faces, sector_faces))
        total_ids = np.vstack((total_ids, sector_ids))
        x_offsets_total = np.vstack((x_offsets_total, x_offsets))
        y_offsets_total = np.vstack((y_offsets_total, y_offsets))
        z_offsets_total = np.vstack((z_offsets_total, z_offsets))
    max_vertex_index = total_vertices.shape[0]
    max_face_index = total_ids.shape[0]

    print('\n Now writing the off file...\n')

    write_to_off_file(
        "DREAM_All.off",
        total_vertices.shape[0],
        total_faces.shape[0],
        total_vertices,
        total_faces,
    )

    print(' Now writing the NeXus file...\n')

    write_to_nexus_file(
        "DREAM_All.nxs",
        det_number,
        total_vertices.shape[0],
        total_vertices,
        total_faces,
        total_ids,
        w_number,
        st_number,
        x_offsets_total,
        y_offsets_total,
        z_offsets_total,
    )

print('Executed successfully!')
