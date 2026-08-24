package org.benchmark.filter4;

import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.io.cif.CifStructureConverter;
import org.biojava.nbio.structure.contact.Contact;
import org.biojava.nbio.structure.contact.Grid;
import org.biojava.nbio.structure.contact.StructureInterface;
import org.biojava.nbio.structure.contact.StructureInterfaceList;
import org.biojava.nbio.structure.xtal.CrystalBuilder;
import org.biojava.nbio.structure.xtal.CrystalTransform;

import javax.vecmath.Matrix4d;
import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.*;
import java.util.*;
import java.util.zip.GZIPInputStream;

public final class Step1Runner {
    private record PairMeta(String pairId,String pdbId,String assemblyId,String modelId,String frameStatus,Set<String> trivial) {}

    public static void main(String[] argv) throws Exception {
        Map<String,String> args = args(argv);
        Path pairFile=Path.of(required(args,"pairs")), atomFile=Path.of(required(args,"atoms"));
        Path mmcifRoot=Path.of(required(args,"mmcif-root")), out=Path.of(required(args,"out-dir"));
        String mode=args.getOrDefault("cell-mode","auto");
        int numCells=Integer.parseInt(args.getOrDefault("num-cells","20"));
        double cutoff=Double.parseDouble(args.getOrDefault("cutoff","6.0"));
        boolean compare=Boolean.parseBoolean(args.getOrDefault("compare-auto-exhaustive","false"));
        Set<String> referencePdbs = args.containsKey("reference-pdbs") ?
                new HashSet<>(Files.readAllLines(Path.of(args.get("reference-pdbs")))) : Set.of();
        Files.createDirectories(out);
        Map<String,PairMeta> meta=readPairs(pairFile);
        Map<String,List<double[]>> ligand=new HashMap<>(), pocket=new HashMap<>(); readAtoms(atomFile,ligand,pocket);
        Map<String,List<CrystalSearchTarget>> byPdb=new TreeMap<>();
        for(PairMeta p:meta.values()) byPdb.computeIfAbsent(p.pdbId,k->new ArrayList<>()).add(new CrystalSearchTarget(
                p.pairId,p.pdbId,p.assemblyId,p.modelId,p.frameStatus,p.trivial,
                ligand.getOrDefault(p.pairId,List.of()),pocket.getOrDefault(p.pairId,List.of())));
        try(BufferedWriter inv=writer(out.resolve("pair_inventory.tsv")); BufferedWriter hits=writer(out.resolve("hits.tsv"));
            BufferedWriter src=writer(out.resolve("source_objects.tsv")); BufferedWriter val=writer(out.resolve("lattice_validation.tsv"))) {
            inv.write("candidate_pair_id\tpdb_id\tstep1_status\ttarget_frame_status\tn_source_objects\tn_raw_lattice_candidates\tn_bbox_pass_candidates\tn_unique_neighbor_instances\tn_ligand_6A_instances\tn_pocket_6A_instances\tn_ligand_and_pocket_6A_instances\tn_trivial_self_instances_skipped\truntime_ms\terror_reason\n");
            hits.write("candidate_pair_id\tpdb_id\tassembly_id\tmodel_id\tcrystal_instance_id\tcrystal_instance_key\tsource_object_id\tsource_object_type\tsource_entity_id\tsource_label_asym_id\tsource_auth_asym_id\tsource_comp_id\tsource_residue_id\tsymmetry_operation_id\tsymmetry_operation_string\tcell_h\tcell_k\tcell_l\tR11\tR12\tR13\tR21\tR22\tR23\tR31\tR32\tR33\ttx\tty\ttz\tmin_distance_to_ligand\tmin_distance_to_pocket\tn_atom_pairs_within_6A_ligand\tn_atom_pairs_within_6A_pocket\ttouches_ligand_6A\ttouches_pocket_6A\tsource_atom_count\tsource_heavy_atom_count\ttarget_frame_status\n");
            src.write("pdb_id\tmodel_id\tsource_object_id\tsource_object_type\tentity_id\tlabel_asym_id\tauth_asym_id\tcomp_id\tresidue_id\tsource_atom_count\tsource_heavy_atom_count\n");
            val.write("pdb_id\tparse_status\tsymmetry_operation_count\ttransform_sample_count\ttransform_max_deviation\tauto_exhaustive_compared_pairs\tauto_exhaustive_missed_instances\tauto_exhaustive_extra_instances\tbiojava_reference_interface_count\tbiojava_reference_missed_placements\tvalidation_status\terror_reason\n");
            for(var entry:byPdb.entrySet()) processPdb(entry.getKey(),entry.getValue(),mmcifRoot,out,mode,numCells,cutoff,compare,referencePdbs.contains(entry.getKey()),inv,hits,src,val);
        }
    }

    private static void processPdb(String pdb,List<CrystalSearchTarget> targets,Path root,Path out,String mode,int numCells,double cutoff,boolean compare,boolean runReference,
                                   BufferedWriter inv,BufferedWriter hitOut,BufferedWriter srcOut,BufferedWriter valOut) throws IOException {
        Path cif=root.resolve(pdb.substring(1,3)).resolve(pdb+".cif.gz");
        try {
            FileParsingParameters params=new FileParsingParameters(); params.setAlignSeqRes(false); params.setCreateAtomBonds(false); params.setParseBioAssembly(false);
            Structure structure;
            try(InputStream in=new GZIPInputStream(Files.newInputStream(cif))) { structure=CifStructureConverter.fromInputStream(in,params); }
            List<CrystalSourceObject> objects=sourceObjects(pdb,structure);
            for(CrystalSourceObject o:objects) srcOut.write(tsv(o.pdbId,o.modelId,o.objectId,o.type,o.entityId,o.labelAsymId,o.authAsymId,o.compId,o.residueId,o.atomCount,o.heavyAtoms.length));
            LigandCrystalBuilder builder=new LigandCrystalBuilder(structure,objects,cutoff,numCells);
            long maxSamples=0; double maxDev=0; long missed=0,extra=0,compared=0;
            for(CrystalSearchTarget target:targets) {
                long start=System.nanoTime();
                LigandCrystalBuilder.SearchResult result=builder.search(target,mode);
                if(compare) {
                    LigandCrystalBuilder.SearchResult a=mode.equals("auto")?result:builder.search(target,"auto");
                    LigandCrystalBuilder.SearchResult e=mode.equals("exhaustive")?result:builder.search(target,"exhaustive");
                    Set<String> ak=keys(a.hits()), ek=keys(e.hits());
                    Set<String> m=new HashSet<>(ek);m.removeAll(ak);missed+=m.size();
                    Set<String> x=new HashSet<>(ak);x.removeAll(ek);extra+=x.size(); compared++;
                }
                int lig=0,poc=0,both=0;
                for(CrystalNeighborHit h:result.hits()) { boolean l=h.ligandContactCount()>0,q=h.pocketContactCount()>0; if(l)lig++;if(q)poc++;if(l&&q)both++; writeHit(hitOut,h); }
                String status=result.hits().isEmpty()?"NO_NEIGHBOR_WITHIN_6A":"SUCCESS";
                long ms=(System.nanoTime()-start)/1_000_000;
                inv.write(tsv(target.pairId,pdb,status,target.targetFrameStatus,objects.size(),result.rawCandidates(),result.bboxPass(),result.hits().size(),lig,poc,both,result.trivialSkipped(),ms,""));
                maxSamples+=result.hits().size(); maxDev=Math.max(maxDev,result.transformMaxDeviation());
            }
            long[] reference = runReference ? compareBioJavaReference(structure,objects,cutoff) : new long[]{0,0};
            String validation=(maxDev<=1e-5 && missed==0 && reference[1]==0)?"PASS":"FAIL";
            valOut.write(tsv(pdb,"PARSE_SUCCESS",builder.symmetryOperationCount(),maxSamples,maxDev,compared,missed,extra,reference[0],reference[1],validation,""));
        } catch(Exception ex) {
            String error=(ex.getClass().getSimpleName()+":"+String.valueOf(ex.getMessage())).replace('\t',' ').replace('\n',' ');
            for(CrystalSearchTarget target:targets) inv.write(tsv(target.pairId,pdb,"STRUCTURE_PARSE_ERROR",target.targetFrameStatus,0,0,0,0,0,0,0,0,0,error));
            valOut.write(tsv(pdb,"PARSE_ERROR",0,0,"",0,0,0,0,0,"FAIL",error));
        }
    }

    private static long[] compareBioJavaReference(Structure structure,List<CrystalSourceObject> objects,double cutoff) {
        StructureInterfaceList refs = new CrystalBuilder(structure).getUniqueInterfaces(cutoff);
        Map<String,Chain> chains = new HashMap<>();
        for(Chain c:structure.getPolyChains(0)){chains.put(clean(c.getName()),c);chains.putIfAbsent(clean(c.getId()),c);}
        Map<String,CrystalSourceObject> sources = new HashMap<>();
        for(CrystalSourceObject o:objects)if(o.type==CrystalSourceObject.Type.POLYMER){sources.put(clean(o.authAsymId),o);sources.putIfAbsent(clean(o.labelAsymId),o);}
        long total=0,missed=0;
        Matrix4d[] ops=structure.getCrystallographicInfo().getTransformationsOrthonormal();
        for(StructureInterface ref:refs){
            total++;
            String iName=ref.getMoleculeIds().getFirst(),jName=ref.getMoleculeIds().getSecond();
            Chain target=chains.get(iName);CrystalSourceObject source=sources.get(jName);CrystalTransform ct=ref.getTransforms().getSecond();
            if(target==null||source==null||ct.getTransformId()<0||ct.getTransformId()>=ops.length){missed++;continue;}
            Matrix4d m=new Matrix4d(ops[ct.getTransformId()]);
            javax.vecmath.Point3i q=ct.getCrystalTranslation();
            javax.vecmath.Vector3d shift=new javax.vecmath.Vector3d(q.x,q.y,q.z);structure.getCrystallographicInfo().getCrystalCell().transfToOrthonormal(shift);
            m.m03+=shift.x;m.m13+=shift.y;m.m23+=shift.z;
            Atom[] targetAtoms=heavyChainAtoms(target),image=transformed(source.heavyAtoms,m);
            if(targetAtoms.length==0||image.length==0){missed++;continue;}
            Grid grid=new Grid(cutoff);grid.addAtoms(targetAtoms,image);List<Contact> contacts=grid.getIndicesContacts();
            if(contacts.isEmpty())missed++;
        }
        return new long[]{total,missed};
    }

    private static Atom[] heavyChainAtoms(Chain c){List<Atom>a=new ArrayList<>();int[]n={0};for(Group g:c.getAtomGroups())collect(g,a,n);return a.toArray(Atom[]::new);}
    private static Atom[] transformed(Atom[] source,Matrix4d m){Atom[]out=new Atom[source.length];for(int i=0;i<source.length;i++){javax.vecmath.Point3d p=new javax.vecmath.Point3d(source[i].getX(),source[i].getY(),source[i].getZ());m.transform(p);AtomImpl a=new AtomImpl();a.setName(source[i].getName());a.setElement(source[i].getElement());a.setX(p.x);a.setY(p.y);a.setZ(p.z);out[i]=a;}return out;}

    private static Set<String> keys(List<CrystalNeighborHit> hits){Set<String>s=new HashSet<>();for(var h:hits)s.add(h.key().toString());return s;}

    private static List<CrystalSourceObject> sourceObjects(String pdb,Structure s) {
        List<CrystalSourceObject> out=new ArrayList<>();
        for(Chain c:s.getPolyChains(0)) {
            String entity=entity(c), label=clean(c.getId()), auth=clean(c.getName()); List<Atom> heavy=new ArrayList<>(); int[] total={0};
            for(Group g:c.getAtomGroups()) collect(g,heavy,total);
            if(!heavy.isEmpty()) {String id="POLYMER|1|"+entity+"|"+label+"|"+auth;out.add(new CrystalSourceObject(pdb,"1",id,CrystalSourceObject.Type.POLYMER,entity,label,auth,"","",total[0],heavy.toArray(Atom[]::new)));}
        }
        for(Chain c:s.getNonPolyChains(0)) {
            String entity=entity(c), label=clean(c.getId()), auth=clean(c.getName());
            for(Group g:c.getAtomGroups()) {
                String comp=clean(g.getPDBName()).toUpperCase(Locale.ROOT); if(comp.equals("HOH")||comp.equals("DOD"))continue;
                List<Atom> heavy=new ArrayList<>();int[] total={0};collect(g,heavy,total);if(heavy.isEmpty())continue;
                ResidueNumber rn=g.getResidueNumber();String authSeq=rn==null||rn.getSeqNum()==null?"":rn.getSeqNum().toString();String ins=rn==null||rn.getInsCode()==null?"":rn.getInsCode().toString();
                String residue=authSeq+"||"+ins+"|"+comp;String id="NONPOLYMER|1|"+entity+"|"+label+"|"+auth+"|"+comp+"|"+authSeq+"||"+ins;
                out.add(new CrystalSourceObject(pdb,"1",id,CrystalSourceObject.Type.NONPOLYMER,entity,label,auth,comp,residue,total[0],heavy.toArray(Atom[]::new)));
            }
        }
        return out;
    }

    private static void collect(Group g,List<Atom> heavy,int[] total){Set<String>seen=new HashSet<>();List<Group>groups=new ArrayList<>();groups.add(g);groups.addAll(g.getAltLocs());for(Group q:groups)for(Atom a:q.getAtoms()){String key=a.getPDBserial()+"|"+a.getName()+"|"+a.getX()+"|"+a.getY()+"|"+a.getZ();if(!seen.add(key))continue;total[0]++;Element e=a.getElement();if(e==Element.H||e==Element.D||e==Element.T)continue;heavy.add(a);}}
    private static String entity(Chain c){return c.getEntityInfo()==null?"":Integer.toString(c.getEntityInfo().getMolId());}
    private static String clean(String x){return x==null?"":x.trim();}

    private static void writeHit(BufferedWriter w,CrystalNeighborHit h)throws IOException{Matrix4d m=h.transform();CrystalInstanceKey k=h.key();CrystalSourceObject s=h.source();w.write(tsv(
            h.target().pairId,h.target().pdbId,h.target().assemblyId,h.target().modelId,k.toString(),k.toString(),s.objectId,s.type,s.entityId,s.labelAsymId,s.authAsymId,s.compId,s.residueId,
            k.symmetryOperationId(),h.symmetryOperationString(),k.h(),k.k(),k.l(),m.m00,m.m01,m.m02,m.m10,m.m11,m.m12,m.m20,m.m21,m.m22,m.m03,m.m13,m.m23,
            finite(h.minDistanceLigand()),finite(h.minDistancePocket()),h.ligandContactCount(),h.pocketContactCount(),h.ligandContactCount()>0,h.pocketContactCount()>0,s.atomCount,s.heavyAtoms.length,h.target().targetFrameStatus));}
    private static Object finite(double x){return Double.isFinite(x)?x:"";}

    private static Map<String,PairMeta> readPairs(Path path)throws IOException{Map<String,PairMeta>out=new LinkedHashMap<>();try(BufferedReader r=Files.newBufferedReader(path)){String[]h=r.readLine().split("\\t",-1);Map<String,Integer>ix=index(h);for(String line;(line=r.readLine())!=null;){String[]v=line.split("\\t",-1);Set<String>tr=new HashSet<>();if(!get(v,ix,"trivial_source_object_ids").isBlank())tr.addAll(Arrays.asList(get(v,ix,"trivial_source_object_ids").split(";")));PairMeta p=new PairMeta(get(v,ix,"candidate_pair_id"),get(v,ix,"pdb_id"),get(v,ix,"assembly_id"),get(v,ix,"model_id"),get(v,ix,"target_frame_status"),tr);out.put(p.pairId,p);}}return out;}
    private static void readAtoms(Path path,Map<String,List<double[]>>lig,Map<String,List<double[]>>poc)throws IOException{try(BufferedReader r=Files.newBufferedReader(path)){String[]h=r.readLine().split("\\t",-1);Map<String,Integer>ix=index(h);for(String line;(line=r.readLine())!=null;){String[]v=line.split("\\t",-1);String id=get(v,ix,"candidate_pair_id"),kind=get(v,ix,"target_kind");double[]xyz={Double.parseDouble(get(v,ix,"x")),Double.parseDouble(get(v,ix,"y")),Double.parseDouble(get(v,ix,"z"))};(kind.equals("LIGAND")?lig:poc).computeIfAbsent(id,k->new ArrayList<>()).add(xyz);}}}
    private static Map<String,Integer>index(String[]h){Map<String,Integer>m=new HashMap<>();for(int i=0;i<h.length;i++)m.put(h[i],i);return m;}
    private static String get(String[]v,Map<String,Integer>i,String k){return v[i.get(k)];}
    private static BufferedWriter writer(Path p)throws IOException{return Files.newBufferedWriter(p,StandardCharsets.UTF_8);}
    private static String tsv(Object...v){StringBuilder b=new StringBuilder();for(int i=0;i<v.length;i++){if(i>0)b.append('\t');String s=String.valueOf(v[i]);b.append(s.replace('\t',' ').replace('\n',' '));}return b.append('\n').toString();}
    private static Map<String,String>args(String[]a){Map<String,String>m=new HashMap<>();for(int i=0;i<a.length;i+=2)m.put(a[i].replaceFirst("^--",""),a[i+1]);return m;}
    private static String required(Map<String,String>m,String k){if(!m.containsKey(k))throw new IllegalArgumentException("missing --"+k);return m.get(k);}
}
