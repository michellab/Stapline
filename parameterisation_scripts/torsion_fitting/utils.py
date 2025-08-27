from typing import Any, Sequence, Tuple, Union, Iterable, Dict, Optional, List
import re

from IPython.display import SVG
from openff.toolkit.topology import FrozenMolecule

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.rdDepictor import Compute2DCoords
from rdkit.Geometry.rdGeometry import Point2D

from openff.fragmenter.fragment import Fragment

from openff.toolkit.topology import Molecule
from openff.toolkit.utils import RDKitToolkitWrapper

__all__ = [
    "DEFAULT_WIDTH",
    "DEFAULT_HEIGHT",
    "draw_molecule",
    "widgetify",
    "display_in_grid",
    "visualize_grid",
    "smiles_to_smarts",
    "depict_fragment",
]

DEFAULT_WIDTH = 300
DEFAULT_HEIGHT = 200

Color = Iterable[float]
BondIndices = Tuple[int, int]


def draw_molecule(
    molecule: Union[Molecule, Chem.Mol],
    image_width: int = DEFAULT_WIDTH,
    image_height: int = DEFAULT_HEIGHT,
    highlight_atoms: Optional[Union[List[int], Dict[int, Color]]] = None,
    highlight_bonds: Optional[
        Union[List[BondIndices], Dict[BondIndices, Color]]
    ] = None,
    atom_notes: Optional[Dict[int, str]] = None,
    bond_notes: Optional[Dict[BondIndices, str]] = None,
    emphasize_atoms: Optional[List[int]] = None,
    explicit_hydrogens: Optional[bool] = None,
    color_by_element: Optional[bool] = None,
) -> SVG:
    """Draw a molecule

    Parameters
    ==========

    molecule
        The molecule to draw.
    image_width
        The width of the resulting image in pixels.
    image_height
        The height of the resulting image in pixels.
    highlight_atoms
        A list of atom indices to highlight, or a map from indices to colors.
        Colors should be given as triplets of floats between 0.0 and 1.0.
    highlight_bonds
        A list of pairs of atom indices indicating bonds to highlight, or a map
        from index pairs to colors. Colors should be given as triplets of floats
        between 0.0 and 1.0.
    atom_notes
        A map from atom indices to a string that should be printed near the
        atom.
    bond_notes
        A map from atom index pairs to a string that should be printed near the
        bond.
    emphasize_atoms
        A list of atom indices to emphasize by drawing other atoms (and their
        bonds) in light grey.
    explicit_hydrogens
        If ``False``, allow uncharged monovalent hydrogens to be hidden. If
        ``True``, make all hydrogens explicit. If ``None``, defer to the
        provided molecule.
    color_by_element
        If True, color heteroatoms according to their element; if False, color
        atoms and bonds monochromatically. By default, uses black and white when
        highlight_atoms or highlight_bonds is provided, and color otherwise.

    Raises
    ======

    KeyError
        When an atom or bond in highlight_atoms or highlight_bonds is missing
        from the image, including when it is present in the molecule but hidden.
    """
    # We're working in RDKit
    if isinstance(molecule, FrozenMolecule):
        rdmol = molecule.to_rdkit()
    else:
        rdmol = Chem.mol(molecule)

    # Process color_by_element argument
    if color_by_element is None:
        color_by_element = highlight_atoms is None and highlight_bonds is None

    if color_by_element:
        set_atom_palette = lambda draw_options: draw_options.useDefaultAtomPalette()
    else:
        set_atom_palette = lambda draw_options: draw_options.useBWAtomPalette()

    # Process explicit_hydrogens argument
    # If we need to remove atoms, create a map from the original indices to the
    # new ones.
    if explicit_hydrogens is None:
        idx_map = {i: i for i in range(rdmol.GetNumAtoms())}
    elif explicit_hydrogens:
        idx_map = {i: i for i in range(rdmol.GetNumAtoms())}
        rdmol = Chem.AddHs(rdmol, explicitOnly=True)
    else:
        idx_map = {
            old: new
            for new, old in enumerate(
                a.GetIdx()
                for a in rdmol.GetAtoms()
                if a.GetAtomicNum() != 1 and a.GetMass() != 1
            )
        }
        rdmol = Chem.RemoveHs(rdmol, updateExplicitCount=True)

    # Process highlight_atoms argument for highlightAtoms and highlightAtomColors
    # highlightAtoms takes a list of atom indices
    # highlightAtomColors takes a mapping from atom indices to colors
    if highlight_atoms is None:
        highlight_atoms = []
        highlight_atom_colors = None
    elif isinstance(highlight_atoms, dict):
        highlight_atom_colors = {
            idx_map[i]: tuple(c) for i, c in highlight_atoms.items() if i in idx_map
        }
        highlight_atoms = list(highlight_atoms.keys())
    else:
        highlight_atoms = [idx_map[i] for i in highlight_atoms if i in idx_map]
        highlight_atom_colors = None

    # Process highlight_bonds argument for highlightBonds and highlightBondColors
    # highlightBonds takes a list of bond indices
    # highlightBondColors takes a mapping from bond indices to colors
    if highlight_bonds is None:
        highlight_bonds = []
        highlight_bond_colors = None
    elif isinstance(highlight_bonds, dict):
        highlight_bond_colors = {
            rdmol.GetBondBetweenAtoms(idx_map[i_a], idx_map[i_b]).GetIdx(): tuple(v)
            for (i_a, i_b), v in highlight_bonds.items()
            if i_a in idx_map and i_b in idx_map
        }

        highlight_bonds = list(highlight_bond_colors.keys())
    else:
        highlight_bonds = [
            rdmol.GetBondBetweenAtoms(idx_map[i_a], idx_map[i_b]).GetIdx()
            for i_a, i_b in highlight_bonds
            if i_a in idx_map and i_b in idx_map
        ]
        highlight_bond_colors = None

    # Process bond_notes argument and place notes in the molecule
    if bond_notes is not None:
        for (i_a, i_b), note in bond_notes.items():
            if i_a not in idx_map or i_b not in idx_map:
                continue
            rdbond = rdmol.GetBondBetweenAtoms(idx_map[i_a], idx_map[i_b])
            rdbond.SetProp("bondNote", note)

    # Process atom_notes argument and place notes in the molecule
    if atom_notes is not None:
        for i, note in atom_notes.items():
            if i not in idx_map:
                continue
            rdatom = rdmol.GetAtomWithIdx(idx_map[i])
            rdatom.SetProp("atomNote", note)

    # Fix kekulization so it is the same for all drawn molecules
    Chem.rdmolops.Kekulize(rdmol, clearAromaticFlags=True)

    # Compute 2D coordinates
    Compute2DCoords(rdmol)

    # Construct the drawing object and get a handle to its options
    drawer = Draw.MolDraw2DSVG(image_width, image_height)
    draw_options = drawer.drawOptions()

    # Specify the scale to fit all atoms
    # This is important for emphasize_atoms
    coords_2d = next(rdmol.GetConformers()).GetPositions()[..., (0, 1)]
    drawer.SetScale(
        image_width,
        image_height,
        Point2D(*(coords_2d.min(axis=0) - 1.0)),
        Point2D(*(coords_2d.max(axis=0) + 1.0)),
    )

    # Set the colors used for each element according to the emphasize_atoms and
    # color_by_element arguments
    if emphasize_atoms:
        draw_options.setAtomPalette(
            {i: (0.8, 0.8, 0.8) for i in range(rdmol.GetNumAtoms())}
        )
    else:
        set_atom_palette(draw_options)

    # Draw the molecule
    # Note that if emphasize_atoms is used, this will be the un-emphasized parts
    # of the molecule
    drawer.DrawMolecule(
        rdmol,
        highlightAtoms=highlight_atoms,
        highlightAtomColors=highlight_atom_colors,
        highlightBonds=highlight_bonds,
        highlightBondColors=highlight_bond_colors,
    )

    # Draw an overlapping molecule for the emphasized atoms
    if emphasize_atoms:
        # Set the atom palette according to the color_by_element argument
        set_atom_palette(draw_options)

        # Create a copy of the molecule that removes atoms that aren't emphasized
        emphasized = Chem.rdchem.RWMol(rdmol)
        emphasized.BeginBatchEdit()
        for i in set(idx_map) - set(emphasize_atoms):
            emphasized.RemoveAtom(idx_map[i])
        emphasized.CommitBatchEdit()

        # Draw the molecule. The scale has been fixed and we're re-using the
        # same coordinates, so this will overlap the background molecule
        drawer.DrawMolecule(emphasized)

    # Finalize the SVG
    drawer.FinishDrawing()

    # Return an SVG object that we can view in notebook
    svg_contents = drawer.GetDrawingText()
    return SVG(svg_contents)


def widgetify(obj, layout: Optional[Dict] = None, **kwargs):
    """Convert an object that can display in a notebook to a widget

    Parameters
    ==========
    layout
        The widget layout. See `https://ipywidgets.readthedocs.io/en/latest/examples/Widget%20Layout.html`_
    **kwargs
        Passed on to ``display()``. See `https://ipython.readthedocs.io/en/stable/api/generated/IPython.display.html#IPython.display.display`_
    """
    import ipywidgets
    from IPython.display import display

    if layout is None:
        layout = {}

    widget = ipywidgets.Output(layout=layout)
    with widget:
        display(obj, **kwargs)
    return widget


def display_in_grid(
    objs: Sequence[Any],
    width: int = DEFAULT_WIDTH,
):
    """Display a sequence of objects in a grid

    Parameters
    ==========
    objs
        Sequence of objects to place in the grid
    width
        width of each object, in pixels
    """
    import ipywidgets

    widgets = [widgetify(obj) for obj in objs]
    width_with_padding = width + width // 20
    return ipywidgets.GridBox(
        widgets,
        layout=ipywidgets.Layout(
            grid_template_columns=f"repeat(auto-fill, {width_with_padding}px)"
        ),
    )


def smiles_to_smarts(smiles: str, keep_mappings: Optional[Sequence[int]] = None):
    """
    Convert a fragment smiles to a SMARTS, optionally keeping only some mappings

    Removes unmapped atoms, as they are added to satisfy a valence model and
    may not be present in the parent molecule. Then removes any maps that are
    not given in keep_mappings, and re-numbers the remaining mappings to be
    consecutive and start at 1.

    Parameters
    ==========
    smiles
        The SMILES to process
    keep_mappings
        A list of mapping indexes to keep. The indices will be remapped in order
        to consecutive values starting at 1. If unspecified or ``None``, all
        existing indices will be kept (and remapped if they are not contiguous).
    """
    if keep_mappings is None:
        keep_mappings = sorted(
            [int(m[0]) for m in re.finditer(r"(?<=:)\d+(?=\])", smiles)]
        )

    no_unmapped = re.sub(r"\[[^:\]]+\]", r"", smiles)
    keep_mappings_pattern = "(" + r"|".join([str(i) for i in keep_mappings]) + r")"
    no_other_mapped = re.sub(
        r":(?!" + keep_mappings_pattern + r"\])\d+(?=\])", r"", no_unmapped
    )
    renumbered = no_other_mapped
    for new, old in enumerate(keep_mappings, start=1):
        renumbered = re.sub(r":" + str(old) + r"\]", r":" + str(new) + r"]", renumbered)

    smiles_out = renumbered
    while "()" in smiles_out:
        smiles_out = renumbered.replace("()", "")
    return smiles_out


def depict_fragment(
    parent: Molecule,
    fragment: Fragment,
    width: int = DEFAULT_WIDTH,
    height: int = DEFAULT_HEIGHT,
):
    """Depict a fragment of the parent molecule"""
    bond_smarts = smiles_to_smarts(fragment.smiles, fragment.bond_indices)
    bond_matches = [match for match in parent.chemical_environment_matches(bond_smarts)]
    bond = bond_matches[0]

    fragment_smarts = smiles_to_smarts(fragment.smiles)
    fragment_matches = [
        match
        for match in parent.chemical_environment_matches(
            fragment_smarts, toolkit_registry=RDKitToolkitWrapper()
        )
    ]
    fragment_indices = fragment_matches[0]

    return draw_molecule(
        molecule=parent,
        highlight_bonds=[bond],
        emphasize_atoms=fragment_indices,
        image_width=width,
        image_height=height,
        explicit_hydrogens=False,
    )
