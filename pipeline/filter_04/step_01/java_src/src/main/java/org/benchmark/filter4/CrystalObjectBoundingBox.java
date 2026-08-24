package org.benchmark.filter4;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.contact.BoundingBox;

import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

public final class CrystalObjectBoundingBox {
    private final BoundingBox[][] boxes;
    private final BoundingBox[] aggregate;

    public CrystalObjectBoundingBox(Matrix4d[] symmetryOps, java.util.List<CrystalSourceObject> objects) {
        boxes = new BoundingBox[symmetryOps.length][objects.size()];
        aggregate = new BoundingBox[symmetryOps.length];
        for (int op = 0; op < symmetryOps.length; op++) {
            for (int j = 0; j < objects.size(); j++) boxes[op][j] = transformBox(objects.get(j).heavyAtoms, symmetryOps[op]);
            aggregate[op] = union(boxes[op]);
        }
    }

    public BoundingBox get(int op, int objectIndex) { return boxes[op][objectIndex]; }

    public BoundingBox aggregateTranslated(int op, Vector3d translation) {
        BoundingBox b = aggregate[op];
        return new BoundingBox(b.xmin + translation.x, b.xmax + translation.x,
                b.ymin + translation.y, b.ymax + translation.y,
                b.zmin + translation.z, b.zmax + translation.z);
    }

    public BoundingBox translated(int op, int objectIndex, Vector3d translation) {
        BoundingBox b = boxes[op][objectIndex];
        return new BoundingBox(b.xmin + translation.x, b.xmax + translation.x,
                b.ymin + translation.y, b.ymax + translation.y,
                b.zmin + translation.z, b.zmax + translation.z);
    }

    private static BoundingBox transformBox(Atom[] atoms, Matrix4d matrix) {
        Point3d[] points = new Point3d[atoms.length];
        for (int i = 0; i < atoms.length; i++) {
            Point3d p = new Point3d(atoms[i].getX(), atoms[i].getY(), atoms[i].getZ());
            matrix.transform(p);
            points[i] = p;
        }
        return new BoundingBox(points);
    }

    private static BoundingBox union(BoundingBox[] values) {
        if (values.length == 0) throw new IllegalArgumentException("Cannot union an empty bounding-box array");
        double xmin = values[0].xmin, xmax = values[0].xmax;
        double ymin = values[0].ymin, ymax = values[0].ymax;
        double zmin = values[0].zmin, zmax = values[0].zmax;
        for (int i = 1; i < values.length; i++) {
            BoundingBox b = values[i];
            xmin = Math.min(xmin, b.xmin); xmax = Math.max(xmax, b.xmax);
            ymin = Math.min(ymin, b.ymin); ymax = Math.max(ymax, b.ymax);
            zmin = Math.min(zmin, b.zmin); zmax = Math.max(zmax, b.zmax);
        }
        return new BoundingBox(xmin, xmax, ymin, ymax, zmin, zmax);
    }
}
