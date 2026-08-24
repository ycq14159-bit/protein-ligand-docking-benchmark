package org.benchmark.filter4;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.contact.BoundingBox;

import javax.vecmath.Point3d;

public final class CrystalSourceObject {
    public enum Type { POLYMER, NONPOLYMER }

    public final String pdbId;
    public final String modelId;
    public final String objectId;
    public final Type type;
    public final String entityId;
    public final String labelAsymId;
    public final String authAsymId;
    public final String compId;
    public final String residueId;
    public final int atomCount;
    public final Atom[] heavyAtoms;
    public final BoundingBox baseBox;

    public CrystalSourceObject(String pdbId, String modelId, String objectId, Type type,
                               String entityId, String labelAsymId, String authAsymId,
                               String compId, String residueId, int atomCount, Atom[] heavyAtoms) {
        this.pdbId = pdbId;
        this.modelId = modelId;
        this.objectId = objectId;
        this.type = type;
        this.entityId = entityId;
        this.labelAsymId = labelAsymId;
        this.authAsymId = authAsymId;
        this.compId = compId;
        this.residueId = residueId;
        this.atomCount = atomCount;
        this.heavyAtoms = heavyAtoms;
        Point3d[] points = new Point3d[heavyAtoms.length];
        for (int i = 0; i < heavyAtoms.length; i++) points[i] = new Point3d(heavyAtoms[i].getX(), heavyAtoms[i].getY(), heavyAtoms[i].getZ());
        this.baseBox = new BoundingBox(points);
    }
}
