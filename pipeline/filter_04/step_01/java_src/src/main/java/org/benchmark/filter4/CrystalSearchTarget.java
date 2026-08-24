package org.benchmark.filter4;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.AtomImpl;
import org.biojava.nbio.structure.Element;
import org.biojava.nbio.structure.contact.BoundingBox;

import javax.vecmath.Point3d;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public final class CrystalSearchTarget {
    public final String pairId;
    public final String pdbId;
    public final String assemblyId;
    public final String modelId;
    public final String targetFrameStatus;
    public final Set<String> trivialSourceObjectIds;
    public final Atom[] ligandAtoms;
    public final Atom[] pocketAtoms;
    public final Atom[] unionAtoms;
    public final BoundingBox unionBox;

    public CrystalSearchTarget(String pairId, String pdbId, String assemblyId, String modelId,
                               String targetFrameStatus, Set<String> trivialSourceObjectIds,
                               List<double[]> ligand, List<double[]> pocket) {
        this.pairId = pairId;
        this.pdbId = pdbId;
        this.assemblyId = assemblyId;
        this.modelId = modelId;
        this.targetFrameStatus = targetFrameStatus;
        this.trivialSourceObjectIds = new HashSet<>(trivialSourceObjectIds);
        this.ligandAtoms = makeAtoms(ligand, "L");
        this.pocketAtoms = makeAtoms(pocket, "P");
        List<Atom> all = new ArrayList<>(ligandAtoms.length + pocketAtoms.length);
        java.util.Collections.addAll(all, ligandAtoms);
        java.util.Collections.addAll(all, pocketAtoms);
        this.unionAtoms = all.toArray(Atom[]::new);
        if (unionAtoms.length == 0) throw new IllegalArgumentException("empty target atom set");
        this.unionBox = new BoundingBox(points(unionAtoms));
    }

    private static Atom[] makeAtoms(List<double[]> coords, String prefix) {
        Atom[] atoms = new Atom[coords.size()];
        for (int i = 0; i < coords.size(); i++) {
            double[] xyz = coords.get(i);
            AtomImpl atom = new AtomImpl();
            atom.setName(prefix + i);
            atom.setElement(Element.C);
            atom.setX(xyz[0]); atom.setY(xyz[1]); atom.setZ(xyz[2]);
            atoms[i] = atom;
        }
        return atoms;
    }

    private static Point3d[] points(Atom[] atoms) {
        Point3d[] points = new Point3d[atoms.length];
        for (int i = 0; i < atoms.length; i++) points[i] = new Point3d(atoms[i].getX(), atoms[i].getY(), atoms[i].getZ());
        return points;
    }
}
