package org.benchmark.filter4;

import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.contact.BoundingBox;
import org.biojava.nbio.structure.contact.Contact;
import org.biojava.nbio.structure.contact.Grid;
import org.biojava.nbio.structure.xtal.CrystalCell;
import org.biojava.nbio.structure.xtal.CrystalTransform;

import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;
import java.util.*;

public final class LigandCrystalBuilder {
    public record SearchResult(List<CrystalNeighborHit> hits, long rawCandidates, long bboxPass,
                               long trivialSkipped, double transformMaxDeviation) {}

    private final Structure structure;
    private final List<CrystalSourceObject> objects;
    private final Matrix4d[] ops;
    private final CrystalCell cell;
    private final CrystalObjectBoundingBox boxes;
    private final double cutoff;
    private final int numCells;

    public LigandCrystalBuilder(Structure structure, List<CrystalSourceObject> objects, double cutoff, int numCells) {
        this.structure = structure;
        this.objects = objects;
        this.cutoff = cutoff;
        this.numCells = numCells;
        PDBCrystallographicInfo info = structure.getCrystallographicInfo();
        if (info == null || info.getCrystalCell() == null || info.getSpaceGroup() == null)
            throw new IllegalArgumentException("BioJava crystal metadata unavailable");
        this.cell = info.getCrystalCell();
        this.ops = info.getTransformationsOrthonormal();
        this.boxes = new CrystalObjectBoundingBox(ops, objects);
    }

    public int symmetryOperationCount() { return ops.length; }

    public SearchResult search(CrystalSearchTarget target, String mode) {
        return mode.equals("exhaustive") ? searchExhaustive(target) : searchAuto(target);
    }

    private SearchResult searchAuto(CrystalSearchTarget target) {
        List<CrystalNeighborHit> hits = new ArrayList<>();
        Set<CrystalInstanceKey> seen = new HashSet<>();
        long raw = 0, bbox = 0, trivial = 0;
        double maxTransformDeviation = 0.0;
        for (int op = 0; op < ops.length; op++) {
            String opString = new CrystalTransform(structure.getCrystallographicInfo().getSpaceGroup(), op).toXYZString();
            for (int objectIndex = 0; objectIndex < objects.size(); objectIndex++) {
                CrystalSourceObject source = objects.get(objectIndex);
                int[] ranges = autoRanges(target.unionBox, boxes.get(op, objectIndex));
                for (int h = ranges[0]; h <= ranges[1]; h++) for (int k = ranges[2]; k <= ranges[3]; k++) for (int l = ranges[4]; l <= ranges[5]; l++) {
                    raw++;
                    CrystalInstanceKey key = new CrystalInstanceKey(source.objectId, op, h, k, l);
                    if (!seen.add(key)) throw new IllegalStateException("duplicate CrystalInstanceKey " + key);
                    if (op == 0 && h == 0 && k == 0 && l == 0 && target.trivialSourceObjectIds.contains(source.objectId)) {
                        trivial++;
                        continue;
                    }
                    Vector3d shift = cellShift(h, k, l);
                    BoundingBox imageBox = boxes.translated(op, objectIndex, shift);
                    if (!target.unionBox.overlaps(imageBox, cutoff)) continue;
                    bbox++;
                    Matrix4d transform = new Matrix4d(ops[op]);
                    transform.m03 += shift.x; transform.m13 += shift.y; transform.m23 += shift.z;
                    Atom[] image = transformAtoms(source.heavyAtoms, transform);
                    ContactSummary ligand = contacts(target.ligandAtoms, image);
                    ContactSummary pocket = contacts(target.pocketAtoms, image);
                    if (ligand.count == 0 && pocket.count == 0) continue;
                    maxTransformDeviation = Math.max(maxTransformDeviation, validateTransform(source.heavyAtoms, image, transform));
                    hits.add(new CrystalNeighborHit(target, source, key, opString, transform,
                            ligand.minimum, pocket.minimum, ligand.count, pocket.count));
                }
            }
        }
        return new SearchResult(hits, raw, bbox, trivial, maxTransformDeviation);
    }

    private SearchResult searchExhaustive(CrystalSearchTarget target) {
        List<CrystalNeighborHit> hits=new ArrayList<>();Set<CrystalInstanceKey>seen=new HashSet<>();
        long raw=(long)objects.size()*ops.length*(2L*numCells+1)*(2L*numCells+1)*(2L*numCells+1),bbox=0,trivial=0;double maxDev=0;
        for(int op=0;op<ops.length;op++){
            String opString=new CrystalTransform(structure.getCrystallographicInfo().getSpaceGroup(),op).toXYZString();
            for(int h=-numCells;h<=numCells;h++)for(int k=-numCells;k<=numCells;k++)for(int l=-numCells;l<=numCells;l++){
                Vector3d shift=cellShift(h,k,l);
                if(!target.unionBox.overlaps(boxes.aggregateTranslated(op,shift),cutoff))continue;
                for(int objectIndex=0;objectIndex<objects.size();objectIndex++){
                    CrystalSourceObject source=objects.get(objectIndex);CrystalInstanceKey key=new CrystalInstanceKey(source.objectId,op,h,k,l);
                    if(!seen.add(key))throw new IllegalStateException("duplicate CrystalInstanceKey "+key);
                    if(op==0&&h==0&&k==0&&l==0&&target.trivialSourceObjectIds.contains(source.objectId)){trivial++;continue;}
                    if(!target.unionBox.overlaps(boxes.translated(op,objectIndex,shift),cutoff))continue;bbox++;
                    Matrix4d transform=new Matrix4d(ops[op]);transform.m03+=shift.x;transform.m13+=shift.y;transform.m23+=shift.z;
                    Atom[] image=transformAtoms(source.heavyAtoms,transform);ContactSummary ligand=contacts(target.ligandAtoms,image),pocket=contacts(target.pocketAtoms,image);
                    if(ligand.count==0&&pocket.count==0)continue;
                    maxDev=Math.max(maxDev,validateTransform(source.heavyAtoms,image,transform));
                    hits.add(new CrystalNeighborHit(target,source,key,opString,transform,ligand.minimum,pocket.minimum,ligand.count,pocket.count));
                }
            }
        }
        return new SearchResult(hits,raw,bbox,trivial,maxDev);
    }

    private int[] autoRanges(BoundingBox target, BoundingBox source) {
        double[] tf = fractionalBounds(expand(target, cutoff));
        double[] sf = fractionalBounds(source);
        int h0 = (int)Math.ceil(tf[0] - sf[1] - 1e-10), h1 = (int)Math.floor(tf[1] - sf[0] + 1e-10);
        int k0 = (int)Math.ceil(tf[2] - sf[3] - 1e-10), k1 = (int)Math.floor(tf[3] - sf[2] + 1e-10);
        int l0 = (int)Math.ceil(tf[4] - sf[5] - 1e-10), l1 = (int)Math.floor(tf[5] - sf[4] + 1e-10);
        return new int[]{h0,h1,k0,k1,l0,l1};
    }

    private double[] fractionalBounds(BoundingBox b) {
        double[] out = {Double.POSITIVE_INFINITY,Double.NEGATIVE_INFINITY,Double.POSITIVE_INFINITY,
                Double.NEGATIVE_INFINITY,Double.POSITIVE_INFINITY,Double.NEGATIVE_INFINITY};
        for (double x : new double[]{b.xmin,b.xmax}) for (double y : new double[]{b.ymin,b.ymax}) for (double z : new double[]{b.zmin,b.zmax}) {
            Point3d p = new Point3d(x,y,z); cell.transfToCrystal(p);
            out[0]=Math.min(out[0],p.x); out[1]=Math.max(out[1],p.x);
            out[2]=Math.min(out[2],p.y); out[3]=Math.max(out[3],p.y);
            out[4]=Math.min(out[4],p.z); out[5]=Math.max(out[5],p.z);
        }
        return out;
    }

    private static BoundingBox expand(BoundingBox b, double c) {
        return new BoundingBox(b.xmin-c,b.xmax+c,b.ymin-c,b.ymax+c,b.zmin-c,b.zmax+c);
    }

    private Vector3d cellShift(int h, int k, int l) {
        Vector3d v = new Vector3d(h,k,l); cell.transfToOrthonormal(v); return v;
    }

    private static Atom[] transformAtoms(Atom[] source, Matrix4d matrix) {
        Atom[] out = new Atom[source.length];
        for (int i=0;i<source.length;i++) {
            Point3d p = new Point3d(source[i].getX(),source[i].getY(),source[i].getZ()); matrix.transform(p);
            AtomImpl a = new AtomImpl(); a.setName(source[i].getName()); a.setElement(source[i].getElement());
            a.setX(p.x); a.setY(p.y); a.setZ(p.z); out[i]=a;
        }
        return out;
    }

    private record ContactSummary(int count, double minimum) {}
    private ContactSummary contacts(Atom[] target, Atom[] image) {
        if (target.length == 0 || image.length == 0) return new ContactSummary(0, Double.NaN);
        Grid grid = new Grid(cutoff); grid.addAtoms(target, image);
        List<Contact> contacts = grid.getIndicesContacts();
        double min = Double.POSITIVE_INFINITY;
        for (Contact contact : contacts) min = Math.min(min, contact.getDistance());
        return new ContactSummary(contacts.size(), contacts.isEmpty()?Double.NaN:min);
    }

    private static double validateTransform(Atom[] source, Atom[] image, Matrix4d m) {
        int n = Math.min(3, source.length); double max=0;
        for (int i=0;i<n;i++) {
            double x=m.m00*source[i].getX()+m.m01*source[i].getY()+m.m02*source[i].getZ()+m.m03;
            double y=m.m10*source[i].getX()+m.m11*source[i].getY()+m.m12*source[i].getZ()+m.m13;
            double z=m.m20*source[i].getX()+m.m21*source[i].getY()+m.m22*source[i].getZ()+m.m23;
            double d=Math.sqrt(Math.pow(x-image[i].getX(),2)+Math.pow(y-image[i].getY(),2)+Math.pow(z-image[i].getZ(),2));
            max=Math.max(max,d);
        }
        return max;
    }
}
