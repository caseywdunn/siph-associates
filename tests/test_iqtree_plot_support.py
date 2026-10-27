"""Check that rooting changes do not move IQ-TREE support between edges."""

import io
import itertools
import sys
import unittest
from pathlib import Path

from Bio import Phylo

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "workflow" / "scripts"))
from plot_mycoplasmatales_figures import node_support, root_on_outgroup, split_key


class IQTreePlotSupportTests(unittest.TestCase):
    def test_support_measures_remain_separate(self):
        tree = Phylo.read(io.StringIO("(A:1,B:2)64.2/71;"), "newick")
        self.assertEqual(node_support(tree.root), (64.2, 71.0))
        for invalid in ("0.999", "nan/100", "100/101"):
            tree.root.name, tree.root.confidence = invalid, None
            with self.assertRaises(ValueError):
                node_support(tree.root)

    def test_rooting_preserves_edge_support_and_pairwise_distances(self):
        tree = Phylo.read(io.StringIO(
            "(((A:1,B:2)81/96:3,C:4)83/97:5,(O1:1,O2:2)87/99:6,O3:7);"
        ), "newick")
        tips = frozenset(tip.name for tip in tree.get_terminals())
        before_support = {
            split_key(frozenset(tip.name for tip in clade.get_terminals()), tips):
                node_support(clade)
            for clade in tree.get_nonterminals() if clade is not tree.root
        }
        distances = {(a, b): tree.distance(a, b)
                     for a, b in itertools.combinations(sorted(tips), 2)}
        outgroup = frozenset({"O1", "O2", "O3"})
        root_on_outgroup(tree, outgroup)
        self.assertEqual(
            {frozenset(tip.name for tip in clade.get_terminals())
             for clade in tree.root.clades},
            {outgroup, tips - outgroup},
        )
        after_support = {
            split_key(frozenset(tip.name for tip in clade.get_terminals()), tips):
                node_support(clade)
            for clade in tree.get_nonterminals() if clade is not tree.root
        }
        self.assertEqual(before_support, after_support)
        for (a, b), distance in distances.items():
            self.assertAlmostEqual(distance, tree.distance(a, b))

    def test_outgroup_must_define_an_existing_split(self):
        tree = Phylo.read(io.StringIO("((A:1,B:2)81/96:3,O1:1,O2:2);"), "newick")
        with self.assertRaisesRegex(ValueError, "does not define a split"):
            root_on_outgroup(tree, frozenset({"A", "O1"}))


if __name__ == "__main__":
    unittest.main()
