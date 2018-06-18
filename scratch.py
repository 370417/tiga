from pymongo.mongo_client import MongoClient
import pymongo
import yaml
from bson import ObjectId
from sage.interfaces.magma import Magma

magma = Magma()

C = MongoClient(port=int(27017))
cap = C.curve_automorphisms.passports

helper_code = '''
/* helper functions for orbits */
/* L is list of IDs */
PrtOrbit:=function(L);
   prtstring:="  - [";
   for i in [1..#L-1] do
      prtstring:=prtstring cat "'" cat L[i] cat "'";
      prtstring:=prtstring cat ",";
   end for;
   prtstring:=prtstring cat "'" cat  L[#L] cat "'" cat "]";
   return prtstring; 
end function;
Qi:=function(k,L);
   left_unchanged:=[L[i] : i in [1..k-1]];
   braid:=[L[k]*L[k+1]*L[k]^-1,L[k]];
   right_unchanged:=[L[i]: i in [k+2..#L] ];
   return left_unchanged cat braid cat right_unchanged;
end function;
BrdOrbit:=function(L,r);
   for ell in L do
      for i in [1..r-1] do
         brd:=Qi(i,ell);
         L join:={brd};
      end for;   
   end for;
   return L;
end function;
Strg:=function(L);
   sg:="[";
   for i in [1..#L-1] do
     sg:= sg cat IntegerToString(L[i]) cat ",";
   end for;
   sg:= sg cat IntegerToString(L[#L]) cat "]";
   return sg;
end function;
OuterAutOrbit:=function(L,G);
   for ell in L do
     for g in G do
        L join:={g(ell)};
     end for;
   end for;
   return L;
end function;
OrbitComputeBraid:=function(Vects,r);
   T:={};
   Orbs:=[];
   genvecs:=Vects;   
   while Vects ne {} do
      rand:=Random(Vects);
      T join:={rand};
      
      tempcount:=0;
      Orb:={rand};
      /* t1:=Cputime(); */
      while (tempcount eq 0) do
         sizeOrb:=#Orb;
         Orb:=BrdOrbit(Orb,r); 
         if #Orb eq sizeOrb  then
            tempcount:=1;  /* stops when outer aut doesn't add more */
            N:=genvecs meet Orb;
            Append(~Orbs,N);
         else
            Vects diff:=Orb;
            if #Vects eq 0 then
               N:=genvecs meet Orb;
               Append(~Orbs,N);
               tempcount:=1;/* stop if all in one orbit or got to last */
            end if;
         end if;    
      end while;   
   end while;  /* while Vects not empty */
   return T, Orbs;
end function;
OrbitComputeAut:=function(Vects,A,r);
   /* T will be list of final orbits, A is outer auts */
   T:={};
   Orbs:=[];
   genvecs:=Vects; 
   while Vects ne {} do
      rand:=Random(Vects);
      T join:={rand};
      tempcount:=0;
      Orb:={rand};
      /* t1:=Cputime(); */
      while (tempcount eq 0) do
         sizeOrb:=#Orb;
         Orb:=OuterAutOrbit(Orb,A);
         Orb:=BrdOrbit(Orb,r);
         if #Orb eq sizeOrb  then
            tempcount:=1;  /* stops when outer aut doesn't add more */
            N:=genvecs meet Orb;
            Append(~Orbs,N);
         else
            Vects diff:=Orb;
            if #Vects eq 0 then
                N:=genvecs meet Orb;
                Append(~Orbs,N);
               tempcount:=1;/* stop if all in one orbit or got to last */
            end if;
         end if;    
      end while;   
   end while;  /* while Vects not empty */
   return T, Orbs;
end function;
'''

top_matter = '''
// The results are stored in a list of records called 'data'
// The results for braid equivalence are stored in a list of records called 'braid_data'
RecFormat:=recformat<group, gp_id, signature, gen_vectors, genus, dimension, r, g0, vector_id>;
data:=[];
braid_data:=[];
// Create group as a permutation group, and generate data which is the same for all entries.
gp_id:={group};
H:=SmallGroup(gp_id[1],gp_id[2]);
n:=#H;
LP:=[];
LG:=[g : g in H];
for i in [1..n] do 
    x:=LG[i]; Tx:=[LG[j]*x : j in [1..n]]; permL:=[];
    for j in [1..n] do
        for k in [1..n] do
            if Tx[j] eq LG[k] then
                permL[j]:=k;
                break;
            end if;
        end for;
    end for;
    Append(~LP,permL);
end for;
G:=PermutationGroup<n|LP>;
signature:={signature};
genus:={genus};
r:=#signature-1;
g0:=signature[1];
dim:=r-3;
S:=Sym(gp_id[1]);
'''

action_code = '''
// Here we add an action to data.
gen_vectors:={gen_vectors};
gen_vectors_as_perm:=[S!perm : perm in gen_vectors];
Append(~data, rec<RecFormat | group:=G,
                              gp_id:=gp_id,
                              signature:=signature,
                              gen_vectors:=gen_vectors_as_perm,
                              genus:=genus,
                              dimension:=dim,
                              r:=r,
                              g0:=g0,
                              vector_id:="{_id}">);
'''

braid_action_code = '''
Append(~braid_data, rec<RecFormat | group:=G,
                                    gp_id:=gp_id,
                                    signature:=signature,
                                    gen_vectors:=gen_vectors_as_perm,
                                    genus:=genus,
                                    dimension:=dim,
                                    r:=r,
                                    g0:=g0,
                                    vector_id:="{_id}">);
'''

orbits_code = '''
prtfile:="orbits_magma2.yml";
/* This next command keeps Magma from printing line breaks */
SetColumns(0);
genvecs:=[[* data[i]`gen_vectors,data[i]`vector_id *] : i in [1..#data]];
braid_genvecs:=[[* braid_data[i]`gen_vectors,braid_data[i]`vector_id *] : i in [1..#braid_data]];
/* MAIN CODE */
if #genvecs gt 0 then
   if #genvecs eq 1 then
      PrintFile(prtfile, "topological:");
      PrintFile(prtfile, "-   "); 
     
   else
      B:=AutomorphismGroup(G);     
      f,A,k:=PermutationRepresentation(B);
      h:=Inverse(f); /* Map from A to B */
      aut:= [h(aL): aL in A | not IsInner(h(aL))];   /* Outer Automorphisms */
      Vects:={g[1] : g in genvecs};
      braid_Vects:={g[i] : g in braid_genvecs}
      BrdRep,BrdOrbs:=OrbitComputeBraid(braid_Vects,#signature-1);
      TopRep,TopOrbs:=OrbitComputeAut(Vects,aut,#signature-1); 
      TopOrbsID:=[];
      for j in [1..#TopOrbs] do
         orb:=TopOrbs[j];
         Append(~TopOrbsID,[* *]);
         for orbt in orb do
            for i in [1..#genvecs] do
               if genvecs[i,1] eq orbt then
                  Append(~TopOrbsID[j],genvecs[i,2]);
                  break i;
               end if;
            end for;  
         end for;
      end for;    
      
      BrdOrbsID:=[];
      for j in [1..#BrdOrbs] do
         orb:=BrdOrbs[j];
         Append(~BrdOrbsID,[*  *]);
         for orbt in orb do
            for i in [1..#braid_genvecs] do
               if braid_genvecs[i,1] eq orbt then
                  Append(~BrdOrbsID[j],braid_genvecs[i,2]);
                  break i;
               end if;
            end for;
         end for;
      end for;
      // PrintFile(prtfile,"topological:");
      print "topological:";
      for j in [1..#TopOrbs]  do
         // PrintFile(prtfile,PrtOrbit(TopOrbsID[j]));
         print PrtOrbit(TopOrbsID[j]);
      end for;
      // PrintFile(prtfile,"braid:");
      print "braid:";
      for j in [1..#BrdOrbs] do
         // PrintFile(prtfile,PrtOrbit(BrdOrbsID[j]));
         print PrtOrbit(BrdOrbsID[j]);
      end for;
  end if; /* whether 1 generating vector or more */
end if;   /* whether any generating vectors */
'''

def update_database(magma_output):
    braid_classes = {}
    top_classes = {}
    for vector_ids in magma_output['braid']:
        representative = min(vector_ids)
        for vector_id in vector_ids:
            braid_classes[vector_id] = representative
    for vector_ids in magma_output['topological']:
        representative = min(vector_ids)
        for vector_id in vector_ids:
            top_classes[vector_id] = representative
    for vector_id in braid_classes:
        braid_class = braid_classes[vector_id]
        top_class = top_classes[vector_id]
        if vector_id == braid_class and vector_id == top_class:
            cap.update_one({'_id': ObjectId(vector_id)}, {'$set': {'braid': braid_class, 'topological': top_class, 'braid_rep': True, 'top_rep': True}})
        elif vector_id != braid_class and vector_id == top_class:
            cap.update_one({'_id': ObjectId(vector_id)}, {'$set': {'braid': braid_class, 'topological': top_class, 'braid_rep': False, 'top_rep': True}})
        elif vector_id == braid_class and vector_id != top_class:
            cap.update_one({'_id': ObjectId(vector_id)}, {'$set': {'braid': braid_class, 'topological': top_class, 'braid_rep': True, 'top_rep': False}})    
        else:
            cap.update_one({'_id': ObjectId(vector_id)}, {'$set': {'braid': braid_class, 'topological': top_class, 'braid_rep': False, 'top_rep': False}})

def run_magma(family):
    if family.count() == 1:
        vector_id = family[0]['_id']
        update_database({
            'braid': [[vector_id]],
            'topological': [[vector_id]]})
        return
    family.sort('total_label', pymongo.DESCENDING)
    big_passports = Set()
    code = helper_code + top_matter.format(**family[0])
    for vector in family:
        compute_braid = False
        if vector['passport_label'] in big_passports:
            compute_braid = True
        elif vector['total_label'][-1] != '1':
            big_passports.add(vector['passport_label'])
            compute_braid = True
        code += action_code.format(**vector)
        if compute_braid:
            code += braid_action_code.format(**vector)
    code += orbits_code
    update_database(yaml.load(magma.eval(code)))

for label in cap.find({'genus': 2}).distinct('label'):
    family = cap.find({'label': label}, {
        'label': 1,
        'passport_label': 1,
        'total_label': 1,
        'gen_vectors': 1,
        'group': 1,
        'signature': 1,
        'genus': 1})
    print family[0]['label']
    run_magma(family)

def add_data():
    cap.insert({u'ndim': 3, u'dim': 3, u'passport_label': u'2.2-1.0.2-2-2-2-2-2.1', u'gen_vectors': [[2, 1], [2, 1], [2, 1], [2, 1], [2, 1], [2, 1]], u'jacobian_decomp': [[2, 1, 2]], u'total_label': u'2.2-1.0.2-2-2-2-2-2.1.1', u'g0': 0, u'cc': [1, 1], u'signature': u'[0,2,2,2,2,2,2]', u'con': u'[2,2,2,2,2,2]', u'eqn': [u'y^2=x(x^{4}+a_{1}x^{3}+a_{2}x^{2}+a_{3}x + 1)'], u'label': u'2.2-1.0.2-2-2-2-2-2', u'cyclic_trigonal': False, u'r': 6, u'hyp_involution': [2, 1], u'genus': 2, u'group': u'[2,1]', u'_id': ObjectId('579ba54c75d28e4971d3ed68'), u'hyperelliptic': True, u'group_order': 2})
    cap.insert({u'ndim': 2, u'dim': 2, u'passport_label': u'2.4-2.0.2-2-2-2-2.1', u'gen_vectors': [[2, 1, 4, 3], [2, 1, 4, 3], [2, 1, 4, 3], [3, 4, 1, 2], [4, 3, 2, 1]], u'jacobian_decomp': [[1, 1, 2], [1, 1, 4]], u'total_label': u'2.4-2.0.2-2-2-2-2.1.1', u'g0': 0, u'cc': [1, 1], u'signature': u'[0,2,2,2,2,2]', u'con': u'[2,2,2,3,4]', u'eqn': [], u'label': u'2.4-2.0.2-2-2-2-2', u'cyclic_trigonal': False, u'r': 5, u'hyp_involution': [2, 1, 4, 3], u'genus': 2, u'group': u'[4,2]', u'_id': ObjectId('579ba54d75d28e4971d3ed69'), u'hyperelliptic': True, u'group_order': 4})
    cap.insert({u'ndim': 2, u'dim': 2, u'passport_label': u'2.4-2.0.2-2-2-2-2.2', u'gen_vectors': [[2, 1, 4, 3], [3, 4, 1, 2], [3, 4, 1, 2], [3, 4, 1, 2], [4, 3, 2, 1]], u'jacobian_decomp': [[1, 1, 3], [1, 1, 4]], u'total_label': u'2.4-2.0.2-2-2-2-2.2.1', u'g0': 0, u'cc': [2, 1], u'signature': u'[0,2,2,2,2,2]', u'con': u'[2,3,3,3,4]', u'eqn': [], u'label': u'2.4-2.0.2-2-2-2-2', u'cyclic_trigonal': False, u'r': 5, u'hyp_involution': [3, 4, 1, 2], u'genus': 2, u'group': u'[4,2]', u'_id': ObjectId('579ba54d75d28e4971d3ed6a'), u'hyperelliptic': True, u'group_order': 4})
    cap.insert({u'ndim': 2, u'dim': 2, u'passport_label': u'2.4-2.0.2-2-2-2-2.3', u'gen_vectors': [[2, 1, 4, 3], [3, 4, 1, 2], [4, 3, 2, 1], [4, 3, 2, 1], [4, 3, 2, 1]], u'jacobian_decomp': [[1, 1, 2], [1, 1, 3]], u'total_label': u'2.4-2.0.2-2-2-2-2.3.1', u'g0': 0, u'cc': [3, 1], u'signature': u'[0,2,2,2,2,2]', u'con': u'[2,3,4,4,4]', u'eqn': [], u'label': u'2.4-2.0.2-2-2-2-2', u'cyclic_trigonal': False, u'r': 5, u'hyp_involution': [4, 3, 2, 1], u'genus': 2, u'group': u'[4,2]', u'_id': ObjectId('579ba54e75d28e4971d3ed6b'), u'hyperelliptic': True, u'group_order': 4})
    cap.insert({u'ndim': 1, u'dim': 1, u'passport_label': u'2.8-3.0.2-2-2-4.1', u'gen_vectors': [[2, 1, 4, 3, 6, 5, 8, 7], [3, 4, 1, 2, 7, 8, 5, 6], [5, 6, 8, 7, 1, 2, 4, 3], [8, 7, 5, 6, 4, 3, 1, 2]], u'jacobian_decomp': [[1, 2, 5]], u'total_label': u'2.8-3.0.2-2-2-4.1.1', u'g0': 0, u'cc': [1, 1], u'signature': u'[0,2,2,2,4]', u'con': u'[2,3,4,5]', u'eqn': [u'y^2=x(x^{4}+a_{1}x^{2}+1)', u'y^2=(x^{2}-1)(x^{4}+a_{1}x^{2}+1)'], u'label': u'2.8-3.0.2-2-2-4', u'cyclic_trigonal': False, u'r': 4, u'hyp_involution': [2, 1, 4, 3, 6, 5, 8, 7], u'genus': 2, u'group': u'[8,3]', u'_id': ObjectId('579ba54e75d28e4971d3ed6c'), u'hyperelliptic': True, u'group_order': 8})
    cap.insert({u'ndim': 0, u'dim': 0, u'passport_label': u'2.10-2.0.2-5-10.1', u'gen_vectors': [[6, 7, 8, 9, 10, 1, 2, 3, 4, 5], [2, 3, 4, 5, 1, 7, 8, 9, 10, 6], [10, 6, 7, 8, 9, 5, 1, 2, 3, 4]], u'jacobian_decomp': [[2, 1, 4]], u'total_label': u'2.10-2.0.2-5-10.1.1', u'g0': 0, u'cc': [1, 1], u'signature': u'[0,2,5,10]', u'con': u'[2,3,10]', u'eqn': [u'y^2=x^{5}+1'], u'label': u'2.10-2.0.2-5-10', u'cyclic_trigonal': False, u'r': 3, u'hyp_involution': [6, 7, 8, 9, 10, 1, 2, 3, 4, 5], u'genus': 2, u'group': u'[10,2]', u'_id': ObjectId('579ba54f75d28e4971d3ed6d'), u'hyperelliptic': True, u'group_order': 10})
    cap.insert({u'ndim': 0, u'dim': 0, u'passport_label': u'2.10-2.0.2-5-10.2', u'gen_vectors': [[6, 7, 8, 9, 10, 1, 2, 3, 4, 5], [3, 4, 5, 1, 2, 8, 9, 10, 6, 7], [9, 10, 6, 7, 8, 4, 5, 1, 2, 3]], u'jacobian_decomp': [[2, 1, 4]], u'total_label': u'2.10-2.0.2-5-10.2.1', u'g0': 0, u'cc': [2, 1], u'signature': u'[0,2,5,10]', u'con': u'[2,4,8]', u'eqn': [u'y^2=x^{5}+1'], u'label': u'2.10-2.0.2-5-10', u'cyclic_trigonal': False, u'r': 3, u'hyp_involution': [6, 7, 8, 9, 10, 1, 2, 3, 4, 5], u'genus': 2, u'group': u'[10,2]', u'_id': ObjectId('579ba54f75d28e4971d3ed6e'), u'hyperelliptic': True, u'group_order': 10})
    cap.insert({u'ndim': 0, u'dim': 0, u'passport_label': u'2.10-2.0.2-5-10.3', u'gen_vectors': [[6, 7, 8, 9, 10, 1, 2, 3, 4, 5], [4, 5, 1, 2, 3, 9, 10, 6, 7, 8], [8, 9, 10, 6, 7, 3, 4, 5, 1, 2]], u'jacobian_decomp': [[2, 1, 4]], u'total_label': u'2.10-2.0.2-5-10.3.1', u'g0': 0, u'cc': [3, 1], u'signature': u'[0,2,5,10]', u'con': u'[2,5,9]', u'eqn': [u'y^2=x^{5}+1'], u'label': u'2.10-2.0.2-5-10', u'cyclic_trigonal': False, u'r': 3, u'hyp_involution': [6, 7, 8, 9, 10, 1, 2, 3, 4, 5], u'genus': 2, u'group': u'[10,2]', u'_id': ObjectId('579ba55075d28e4971d3ed6f'), u'hyperelliptic': True, u'group_order': 10})
    cap.insert({u'ndim': 0, u'dim': 0, u'passport_label': u'2.10-2.0.2-5-10.4', u'gen_vectors': [[6, 7, 8, 9, 10, 1, 2, 3, 4, 5], [5, 1, 2, 3, 4, 10, 6, 7, 8, 9], [7, 8, 9, 10, 6, 2, 3, 4, 5, 1]], u'jacobian_decomp': [[2, 1, 4]], u'total_label': u'2.10-2.0.2-5-10.4.1', u'g0': 0, u'cc': [4, 1], u'signature': u'[0,2,5,10]', u'con': u'[2,6,7]', u'eqn': [u'y^2=x^{5}+1'], u'label': u'2.10-2.0.2-5-10', u'cyclic_trigonal': False, u'r': 3, u'hyp_involution': [6, 7, 8, 9, 10, 1, 2, 3, 4, 5], u'genus': 2, u'group': u'[10,2]', u'_id': ObjectId('579ba55075d28e4971d3ed70'), u'hyperelliptic': True, u'group_order': 10})
    cap.insert({u'ndim': 1, u'dim': 1, u'passport_label': u'2.12-4.0.2-2-2-3.1', u'gen_vectors': [[4, 5, 6, 1, 2, 3, 10, 11, 12, 7, 8, 9], [7, 9, 8, 10, 12, 11, 1, 3, 2, 4, 6, 5], [12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1], [2, 3, 1, 5, 6, 4, 8, 9, 7, 11, 12, 10]], u'jacobian_decomp': [[1, 2, 6]], u'total_label': u'2.12-4.0.2-2-2-3.1.1', u'g0': 0, u'cc': [1, 1], u'signature': u'[0,2,2,2,3]', u'con': u'[2,3,4,5]', u'eqn': [u'y^2=x^{6}+a_1x^{3}+1'], u'label': u'2.12-4.0.2-2-2-3', u'cyclic_trigonal': False, u'r': 4, u'hyp_involution': [4, 5, 6, 1, 2, 3, 10, 11, 12, 7, 8, 9], u'genus': 2, u'group': u'[12,4]', u'_id': ObjectId('579ba55175d28e4971d3ed71'), u'hyperelliptic': True, u'group_order': 12})
    cap.insert({u'ndim': 0, u'dim': 0, u'passport_label': u'2.24-8.0.2-4-6.1', u'gen_vectors': [[13, 15, 14, 16, 18, 17, 22, 24, 23, 19, 21, 20, 1, 3, 2, 4, 6, 5, 10, 12, 11, 7, 9, 8], [21, 20, 19, 24, 23, 22, 18, 17, 16, 15, 14, 13, 9, 8, 7, 12, 11, 10, 6, 5, 4, 3, 2, 1], [8, 9, 7, 11, 12, 10, 2, 3, 1, 5, 6, 4, 20, 21, 19, 23, 24, 22, 14, 15, 13, 17, 18, 16]], u'jacobian_decomp': [[1, 2, 8]], u'total_label': u'2.24-8.0.2-4-6.1.1', u'g0': 0, u'cc': [1, 1], u'signature': u'[0,2,4,6]', u'con': u'[4,6,8]', u'eqn': [u'y^2=x^{6}-1'], u'label': u'2.24-8.0.2-4-6', u'cyclic_trigonal': False, u'r': 3, u'hyp_involution': [4, 5, 6, 1, 2, 3, 10, 11, 12, 7, 8, 9, 16, 17, 18, 13, 14, 15, 22, 23, 24, 19, 20, 21], u'genus': 2, u'group': u'[24,8]', u'_id': ObjectId('579ba55175d28e4971d3ed72'), u'hyperelliptic': True, u'group_order': 24})
    cap.insert({u'ndim': 0, u'dim': 0, u'passport_label': u'2.24-8.0.2-4-6.2', u'gen_vectors': [[13, 15, 14, 16, 18, 17, 22, 24, 23, 19, 21, 20, 1, 3, 2, 4, 6, 5, 10, 12, 11, 7, 9, 8], [20, 19, 21, 23, 22, 24, 17, 16, 18, 14, 13, 15, 8, 7, 9, 11, 10, 12, 5, 4, 6, 2, 1, 3], [9, 7, 8, 12, 10, 11, 3, 1, 2, 6, 4, 5, 21, 19, 20, 24, 22, 23, 15, 13, 14, 18, 16, 17]], u'jacobian_decomp': [[1, 2, 8]], u'total_label': u'2.24-8.0.2-4-6.2.1', u'g0': 0, u'cc': [2, 1], u'signature': u'[0,2,4,6]', u'con': u'[4,6,9]', u'eqn': [u'y^2=x^{6}-1'], u'label': u'2.24-8.0.2-4-6', u'cyclic_trigonal': False, u'r': 3, u'hyp_involution': [4, 5, 6, 1, 2, 3, 10, 11, 12, 7, 8, 9, 16, 17, 18, 13, 14, 15, 22, 23, 24, 19, 20, 21], u'genus': 2, u'group': u'[24,8]', u'_id': ObjectId('579ba55275d28e4971d3ed73'), u'hyperelliptic': True, u'group_order': 24})
    cap.insert({u'ndim': 0, u'dim': 0, u'passport_label': u'2.48-29.0.2-3-8.1', u'gen_vectors': [[25, 26, 29, 30, 27, 28, 32, 31, 41, 42, 45, 46, 43, 44, 48, 47, 33, 34, 37, 38, 35, 36, 40, 39, 1, 2, 5, 6, 3, 4, 8, 7, 17, 18, 21, 22, 19, 20, 24, 23, 9, 10, 13, 14, 11, 12, 16, 15], [13, 14, 11, 12, 15, 16, 9, 10, 21, 22, 19, 20, 23, 24, 17, 18, 5, 6, 3, 4, 7, 8, 1, 2, 37, 38, 35, 36, 39, 40, 33, 34, 45, 46, 43, 44, 47, 48, 41, 42, 29, 30, 27, 28, 31, 32, 25, 26], [40, 39, 37, 38, 33, 34, 35, 36, 32, 31, 29, 30, 25, 26, 27, 28, 48, 47, 45, 46, 41, 42, 43, 44, 16, 15, 13, 14, 9, 10, 11, 12, 8, 7, 5, 6, 1, 2, 3, 4, 24, 23, 21, 22, 17, 18, 19, 20]], u'jacobian_decomp': [[1, 2, 4]], u'total_label': u'2.48-29.0.2-3-8.1.1', u'g0': 0, u'cc': [1, 1], u'signature': u'[0,2,3,8]', u'con': u'[3,4,7]', u'eqn': [u'y^2=x(x^4-1)'], u'label': u'2.48-29.0.2-3-8', u'cyclic_trigonal': False, u'r': 3, u'hyp_involution': [2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13, 16, 15, 18, 17, 20, 19, 22, 21, 24, 23, 26, 25, 28, 27, 30, 29, 32, 31, 34, 33, 36, 35, 38, 37, 40, 39, 42, 41, 44, 43, 46, 45, 48, 47], u'genus': 2, u'group': u'[48,29]', u'_id': ObjectId('579ba55275d28e4971d3ed74'), u'hyperelliptic': True, u'group_order': 48})
    cap.insert({u'ndim': 0, u'dim': 0, u'passport_label': u'2.48-29.0.2-3-8.2', u'gen_vectors': [[25, 26, 29, 30, 27, 28, 32, 31, 41, 42, 45, 46, 43, 44, 48, 47, 33, 34, 37, 38, 35, 36, 40, 39, 1, 2, 5, 6, 3, 4, 8, 7, 17, 18, 21, 22, 19, 20, 24, 23, 9, 10, 13, 14, 11, 12, 16, 15], [16, 15, 9, 10, 13, 14, 12, 11, 24, 23, 17, 18, 21, 22, 20, 19, 8, 7, 1, 2, 5, 6, 4, 3, 40, 39, 33, 34, 37, 38, 36, 35, 48, 47, 41, 42, 45, 46, 44, 43, 32, 31, 25, 26, 29, 30, 28, 27], [37, 38, 39, 40, 35, 36, 34, 33, 29, 30, 31, 32, 27, 28, 26, 25, 45, 46, 47, 48, 43, 44, 42, 41, 13, 14, 15, 16, 11, 12, 10, 9, 5, 6, 7, 8, 3, 4, 2, 1, 21, 22, 23, 24, 19, 20, 18, 17]], u'jacobian_decomp': [[1, 2, 4]], u'total_label': u'2.48-29.0.2-3-8.2.1', u'g0': 0, u'cc': [2, 1], u'signature': u'[0,2,3,8]', u'con': u'[3,4,8]', u'eqn': [u'y^2=x(x^4-1)'], u'label': u'2.48-29.0.2-3-8', u'cyclic_trigonal': False, u'r': 3, u'hyp_involution': [2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13, 16, 15, 18, 17, 20, 19, 22, 21, 24, 23, 26, 25, 28, 27, 30, 29, 32, 31, 34, 33, 36, 35, 38, 37, 40, 39, 42, 41, 44, 43, 46, 45, 48, 47], u'genus': 2, u'group': u'[48,29]', u'_id': ObjectId('579ba55375d28e4971d3ed75'), u'hyperelliptic': True, u'group_order': 48})
    cap.insert({u'ndim': 1, u'dim': 1, u'full_auto': u'[12,4]', u'passport_label': u'2.6-1.0.2-2-3-3.1', u'gen_vectors': [[4, 6, 5, 1, 3, 2], [4, 6, 5, 1, 3, 2], [2, 3, 1, 5, 6, 4], [3, 1, 2, 6, 4, 5]], u'jacobian_decomp': [[1, 2, 3]], u'total_label': u'2.6-1.0.2-2-3-3.1.1', u'g0': 0, u'cc': [1, 1], u'eqn': [], u'signH': u'[0,2,2,2,3]', u'con': u'[2,2,3,3]', u'label': u'2.6-1.0.2-2-3-3', u'r': 4, u'genus': 2, u'full_label': u'2.12-4.0.2-2-2-3', u'signature': u'[0,2,2,3,3]', u'group': u'[6,1]', u'_id': ObjectId('579ba55775d28e4971d3ed7d'), u'group_order': 6})
    cap.insert({u'ndim': 1, u'dim': 1, u'full_auto': u'[12,4]', u'passport_label': u'2.6-1.0.2-2-3-3.1', u'gen_vectors': [[4, 6, 5, 1, 3, 2], [6, 5, 4, 3, 2, 1], [3, 1, 2, 6, 4, 5], [3, 1, 2, 6, 4, 5]], u'jacobian_decomp': [[1, 2, 3]], u'total_label': u'2.6-1.0.2-2-3-3.1.2', u'g0': 0, u'cc': [1, 2], u'eqn': [], u'signH': u'[0,2,2,2,3]', u'con': u'[2,2,3,3]', u'label': u'2.6-1.0.2-2-3-3', u'r': 4, u'genus': 2, u'full_label': u'2.12-4.0.2-2-2-3', u'signature': u'[0,2,2,3,3]', u'group': u'[6,1]', u'_id': ObjectId('579ba55775d28e4971d3ed7e'), u'group_order': 6})
    cap.insert({u'ndim': 1, u'dim': 1, u'full_auto': u'[12,4]', u'passport_label': u'2.6-2.0.2-2-3-3.1', u'gen_vectors': [[4, 5, 6, 1, 2, 3], [4, 5, 6, 1, 2, 3], [2, 3, 1, 5, 6, 4], [3, 1, 2, 6, 4, 5]], u'jacobian_decomp': [[2, 1, 4]], u'total_label': u'2.6-2.0.2-2-3-3.1.1', u'g0': 0, u'cc': [1, 1], u'eqn': [], u'signH': u'[0,2,2,2,3]', u'con': u'[2,2,3,4]', u'label': u'2.6-2.0.2-2-3-3', u'r': 4, u'genus': 2, u'full_label': u'2.12-4.0.2-2-2-3', u'signature': u'[0,2,2,3,3]', u'group': u'[6,2]', u'_id': ObjectId('579ba55875d28e4971d3ed7f'), u'group_order': 6})
    cap.insert({u'ndim': 0, u'dim': 0, u'full_auto': u'[48,29]', u'passport_label': u'2.8-4.0.4-4-4.1', u'gen_vectors': [[3, 4, 2, 1, 7, 8, 6, 5], [5, 6, 8, 7, 2, 1, 3, 4], [7, 8, 5, 6, 4, 3, 2, 1]], u'jacobian_decomp': [[2, 1, 5]], u'total_label': u'2.8-4.0.4-4-4.1.1', u'g0': 0, u'cc': [1, 1], u'eqn': [], u'signH': u'[0,2,3,8]', u'con': u'[3,4,5]', u'label': u'2.8-4.0.4-4-4', u'r': 3, u'genus': 2, u'full_label': u'2.48-29.0.2-3-8', u'signature': u'[0,4,4,4]', u'group': u'[8,4]', u'_id': ObjectId('579ba55875d28e4971d3ed80'), u'group_order': 8})
    cap.insert({u'ndim': 0, u'dim': 0, u'full_auto': u'[48,29]', u'passport_label': u'2.8-1.0.2-8-8.1', u'gen_vectors': [[2, 1, 4, 3, 6, 5, 8, 7], [5, 6, 7, 8, 3, 4, 2, 1], [7, 8, 6, 5, 2, 1, 4, 3]], u'jacobian_decomp': [[2, 1, 2]], u'total_label': u'2.8-1.0.2-8-8.1.1', u'g0': 0, u'cc': [1, 1], u'eqn': [], u'signH': u'[0,2,3,8]', u'con': u'[2,5,6]', u'label': u'2.8-1.0.2-8-8', u'r': 3, u'genus': 2, u'full_label': u'2.48-29.0.2-3-8', u'signature': u'[0,2,8,8]', u'group': u'[8,1]', u'_id': ObjectId('579ba55975d28e4971d3ed81'), u'group_order': 8})
    cap.insert({u'ndim': 0, u'dim': 0, u'full_auto': u'[48,29]', u'passport_label': u'2.8-1.0.2-8-8.2', u'gen_vectors': [[2, 1, 4, 3, 6, 5, 8, 7], [6, 5, 8, 7, 4, 3, 1, 2], [8, 7, 5, 6, 1, 2, 3, 4]], u'jacobian_decomp': [[2, 1, 2]], u'total_label': u'2.8-1.0.2-8-8.2.1', u'g0': 0, u'cc': [2, 1], u'eqn': [], u'signH': u'[0,2,3,8]', u'con': u'[2,7,8]', u'label': u'2.8-1.0.2-8-8', u'r': 3, u'genus': 2, u'full_label': u'2.48-29.0.2-3-8', u'signature': u'[0,2,8,8]', u'group': u'[8,1]', u'_id': ObjectId('579ba55975d28e4971d3ed82'), u'group_order': 8})
    cap.insert({u'ndim': 0, u'dim': 0, u'full_auto': u'[24,8]', u'passport_label': u'2.12-1.0.3-4-4.1', u'gen_vectors': [[2, 3, 1, 5, 6, 4, 8, 9, 7, 11, 12, 10], [7, 9, 8, 10, 12, 11, 4, 6, 5, 1, 3, 2], [12, 11, 10, 9, 8, 7, 3, 2, 1, 6, 5, 4]], u'jacobian_decomp': [[2, 1, 5]], u'total_label': u'2.12-1.0.3-4-4.1.1', u'g0': 0, u'cc': [1, 1], u'eqn': [], u'signH': u'[0,2,4,6]', u'con': u'[3,4,5]', u'label': u'2.12-1.0.3-4-4', u'r': 3, u'genus': 2, u'full_label': u'2.24-8.0.2-4-6', u'signature': u'[0,3,4,4]', u'group': u'[12,1]', u'_id': ObjectId('579ba55975d28e4971d3ed83'), u'group_order': 12})
    cap.insert({u'ndim': 0, u'dim': 0, u'full_auto': u'[24,8]', u'passport_label': u'2.12-5.0.2-6-6.1', u'gen_vectors': [[4, 5, 6, 1, 2, 3, 10, 11, 12, 7, 8, 9], [8, 9, 7, 11, 12, 10, 2, 3, 1, 5, 6, 4], [12, 10, 11, 9, 7, 8, 6, 4, 5, 3, 1, 2]], u'jacobian_decomp': [[1, 1, 6], [1, 1, 8]], u'total_label': u'2.12-5.0.2-6-6.1.1', u'g0': 0, u'cc': [1, 1], u'eqn': [], u'signH': u'[0,2,4,6]', u'con': u'[2,9,12]', u'label': u'2.12-5.0.2-6-6', u'r': 3, u'genus': 2, u'full_label': u'2.24-8.0.2-4-6', u'signature': u'[0,2,6,6]', u'group': u'[12,5]', u'_id': ObjectId('579ba55a75d28e4971d3ed84'), u'group_order': 12})
    cap.insert({u'ndim': 0, u'dim': 0, u'full_auto': u'[24,8]', u'passport_label': u'2.12-5.0.2-6-6.2', u'gen_vectors': [[4, 5, 6, 1, 2, 3, 10, 11, 12, 7, 8, 9], [9, 7, 8, 12, 10, 11, 3, 1, 2, 6, 4, 5], [11, 12, 10, 8, 9, 7, 5, 6, 4, 2, 3, 1]], u'jacobian_decomp': [[1, 1, 6], [1, 1, 8]], u'total_label': u'2.12-5.0.2-6-6.2.1', u'g0': 0, u'cc': [2, 1], u'eqn': [], u'signH': u'[0,2,4,6]', u'con': u'[2,10,11]', u'label': u'2.12-5.0.2-6-6', u'r': 3, u'genus': 2, u'full_label': u'2.24-8.0.2-4-6', u'signature': u'[0,2,6,6]', u'group': u'[12,5]', u'_id': ObjectId('579ba55a75d28e4971d3ed85'), u'group_order': 12})
    cap.insert({u'ndim': 0, u'dim': 0, u'full_auto': u'[24,8]', u'passport_label': u'2.12-5.0.2-6-6.3', u'gen_vectors': [[7, 8, 9, 10, 11, 12, 1, 2, 3, 4, 5, 6], [5, 6, 4, 2, 3, 1, 11, 12, 10, 8, 9, 7], [12, 10, 11, 9, 7, 8, 6, 4, 5, 3, 1, 2]], u'jacobian_decomp': [[1, 1, 6], [1, 1, 8]], u'total_label': u'2.12-5.0.2-6-6.3.1', u'g0': 0, u'cc': [3, 1], u'eqn': [], u'signH': u'[0,2,4,6]', u'con': u'[3,7,12]', u'label': u'2.12-5.0.2-6-6', u'r': 3, u'genus': 2, u'full_label': u'2.24-8.0.2-4-6', u'signature': u'[0,2,6,6]', u'group': u'[12,5]', u'_id': ObjectId('579ba55b75d28e4971d3ed86'), u'group_order': 12})
    cap.insert({u'ndim': 0, u'dim': 0, u'full_auto': u'[24,8]', u'passport_label': u'2.12-5.0.2-6-6.4', u'gen_vectors': [[7, 8, 9, 10, 11, 12, 1, 2, 3, 4, 5, 6], [6, 4, 5, 3, 1, 2, 12, 10, 11, 9, 7, 8], [11, 12, 10, 8, 9, 7, 5, 6, 4, 2, 3, 1]], u'jacobian_decomp': [[1, 1, 6], [1, 1, 8]], u'total_label': u'2.12-5.0.2-6-6.4.1', u'g0': 0, u'cc': [4, 1], u'eqn': [], u'signH': u'[0,2,4,6]', u'con': u'[3,8,11]', u'label': u'2.12-5.0.2-6-6', u'r': 3, u'genus': 2, u'full_label': u'2.24-8.0.2-4-6', u'signature': u'[0,2,6,6]', u'group': u'[12,5]', u'_id': ObjectId('579ba55b75d28e4971d3ed87'), u'group_order': 12})
    cap.insert({u'ndim': 0, u'dim': 0, u'full_auto': u'[24,8]', u'passport_label': u'2.12-5.0.2-6-6.5', u'gen_vectors': [[10, 11, 12, 7, 8, 9, 4, 5, 6, 1, 2, 3], [5, 6, 4, 2, 3, 1, 11, 12, 10, 8, 9, 7], [9, 7, 8, 12, 10, 11, 3, 1, 2, 6, 4, 5]], u'jacobian_decomp': [[1, 1, 6], [1, 1, 8]], u'total_label': u'2.12-5.0.2-6-6.5.1', u'g0': 0, u'cc': [5, 1], u'eqn': [], u'signH': u'[0,2,4,6]', u'con': u'[4,7,10]', u'label': u'2.12-5.0.2-6-6', u'r': 3, u'genus': 2, u'full_label': u'2.24-8.0.2-4-6', u'signature': u'[0,2,6,6]', u'group': u'[12,5]', u'_id': ObjectId('579ba55c75d28e4971d3ed88'), u'group_order': 12})
    cap.insert({u'ndim': 0, u'dim': 0, u'full_auto': u'[24,8]', u'passport_label': u'2.12-5.0.2-6-6.6', u'gen_vectors': [[10, 11, 12, 7, 8, 9, 4, 5, 6, 1, 2, 3], [6, 4, 5, 3, 1, 2, 12, 10, 11, 9, 7, 8], [8, 9, 7, 11, 12, 10, 2, 3, 1, 5, 6, 4]], u'jacobian_decomp': [[1, 1, 6], [1, 1, 8]], u'total_label': u'2.12-5.0.2-6-6.6.1', u'g0': 0, u'cc': [6, 1], u'eqn': [], u'signH': u'[0,2,4,6]', u'con': u'[4,8,9]', u'label': u'2.12-5.0.2-6-6', u'r': 3, u'genus': 2, u'full_label': u'2.24-8.0.2-4-6', u'signature': u'[0,2,6,6]', u'group': u'[12,5]', u'_id': ObjectId('579ba55c75d28e4971d3ed89'), u'group_order': 12})
    cap.insert({u'ndim': 0, u'dim': 0, u'full_auto': u'[48,29]', u'passport_label': u'2.24-3.0.3-3-4.1', u'gen_vectors': [[9, 10, 15, 16, 11, 12, 13, 14, 17, 18, 23, 24, 19, 20, 21, 22, 1, 2, 7, 8, 3, 4, 5, 6], [24, 23, 19, 20, 17, 18, 22, 21, 8, 7, 3, 4, 1, 2, 6, 5, 16, 15, 11, 12, 9, 10, 14, 13], [7, 8, 5, 6, 4, 3, 2, 1, 15, 16, 13, 14, 12, 11, 10, 9, 23, 24, 21, 22, 20, 19, 18, 17]], u'jacobian_decomp': [[2, 1, 4]], u'total_label': u'2.24-3.0.3-3-4.1.1', u'g0': 0, u'cc': [1, 1], u'eqn': [], u'signH': u'[0,2,3,8]', u'con': u'[3,4,5]', u'label': u'2.24-3.0.3-3-4', u'r': 3, u'genus': 2, u'full_label': u'2.48-29.0.2-3-8', u'signature': u'[0,3,3,4]', u'group': u'[24,3]', u'_id': ObjectId('579ba55d75d28e4971d3ed8a'), u'group_order': 24})
    cap.insert({u'ndim': 0, u'dim': 0, u'full_auto': u'[48,29]', u'passport_label': u'2.16-8.0.2-4-8.1', u'gen_vectors': [[5, 6, 8, 7, 1, 2, 4, 3, 13, 14, 16, 15, 9, 10, 12, 11], [10, 9, 11, 12, 16, 15, 14, 13, 1, 2, 4, 3, 7, 8, 5, 6], [13, 14, 15, 16, 12, 11, 9, 10, 6, 5, 8, 7, 3, 4, 2, 1]], u'jacobian_decomp': [[1, 2, 6]], u'total_label': u'2.16-8.0.2-4-8.1.1', u'g0': 0, u'cc': [1, 1], u'eqn': [], u'signH': u'[0,2,3,8]', u'con': u'[3,5,6]', u'label': u'2.16-8.0.2-4-8', u'r': 3, u'genus': 2, u'full_label': u'2.48-29.0.2-3-8', u'signature': u'[0,2,4,8]', u'group': u'[16,8]', u'_id': ObjectId('579ba55d75d28e4971d3ed8c'), u'group_order': 16})
    cap.insert({u'ndim': 0, u'dim': 0, u'full_auto': u'[48,29]', u'passport_label': u'2.16-8.0.2-4-8.2', u'gen_vectors': [[5, 6, 8, 7, 1, 2, 4, 3, 13, 14, 16, 15, 9, 10, 12, 11], [9, 10, 12, 11, 15, 16, 13, 14, 2, 1, 3, 4, 8, 7, 6, 5], [14, 13, 16, 15, 11, 12, 10, 9, 5, 6, 7, 8, 4, 3, 1, 2]], u'jacobian_decomp': [[1, 2, 6]], u'total_label': u'2.16-8.0.2-4-8.2.1', u'g0': 0, u'cc': [2, 1], u'eqn': [], u'signH': u'[0,2,3,8]', u'con': u'[3,5,7]', u'label': u'2.16-8.0.2-4-8', u'r': 3, u'genus': 2, u'full_label': u'2.48-29.0.2-3-8', u'signature': u'[0,2,4,8]', u'group': u'[16,8]', u'_id': ObjectId('579ba55e75d28e4971d3ed8d'), u'group_order': 16})
    cap.insert({u'ndim': 0, u'dim': 0, u'full_auto': u'[24,8]', u'passport_label': u'2.6-2.0.3-6-6.2', u'gen_vectors': [[3, 1, 2, 6, 4, 5], [6, 4, 5, 3, 1, 2], [6, 4, 5, 3, 1, 2]], u'jacobian_decomp': [[1, 1, 3], [1, 1, 4]], u'total_label': u'2.6-2.0.3-6-6.2.1', u'g0': 0, u'cc': [2, 1], u'eqn': [], u'signH': u'[0,2,4,6]', u'con': u'[4,6,6]', u'label': u'2.6-2.0.3-6-6', u'r': 3, u'genus': 2, u'full_label': u'2.24-8.0.2-4-6', u'signature': u'[0,3,6,6]', u'group': u'[6,2]', u'_id': ObjectId('579ba55675d28e4971d3ed7c'), u'group_order': 6})
    cap.insert({u'ndim': 1, u'dim': 1, u'full_auto': u'[12,4]', u'passport_label': u'2.3-1.0.3-3-3-3.1', u'gen_vectors': [[2, 3, 1], [2, 3, 1], [3, 1, 2], [3, 1, 2]], u'jacobian_decomp': [[2, 1, 2]], u'total_label': u'2.3-1.0.3-3-3-3.1.1', u'g0': 0, u'cc': [1, 1], u'eqn': [], u'signH': u'[0,2,2,2,3]', u'con': u'[2,2,3,3]', u'label': u'2.3-1.0.3-3-3-3', u'r': 4, u'genus': 2, u'full_label': u'2.12-4.0.2-2-2-3', u'signature': u'[0,3,3,3,3]', u'group': u'[3,1]', u'_id': ObjectId('579ba55d75d28e4971d3ed8b'), u'group_order': 3})
    cap.insert({u'ndim': 0, u'dim': 0, u'full_auto': u'[10,2]', u'passport_label': u'2.5-1.0.5-5-5.4', u'gen_vectors': [[4, 5, 1, 2, 3], [4, 5, 1, 2, 3], [5, 1, 2, 3, 4]], u'jacobian_decomp': [[2, 1, 2]], u'total_label': u'2.5-1.0.5-5-5.4.1', u'g0': 0, u'cc': [4, 1], u'eqn': [], u'signH': u'[0,2,5,10]', u'con': u'[4,4,5]', u'label': u'2.5-1.0.5-5-5', u'r': 3, u'genus': 2, u'full_label': u'2.10-2.0.2-5-10', u'signature': u'[0,5,5,5]', u'group': u'[5,1]', u'_id': ObjectId('579ba55575d28e4971d3ed7a'), u'group_order': 5})
    cap.insert({u'ndim': 0, u'dim': 0, u'full_auto': u'[24,8]', u'passport_label': u'2.6-2.0.3-6-6.1', u'gen_vectors': [[2, 3, 1, 5, 6, 4], [5, 6, 4, 2, 3, 1], [5, 6, 4, 2, 3, 1]], u'jacobian_decomp': [[1, 1, 3], [1, 1, 4]], u'total_label': u'2.6-2.0.3-6-6.1.1', u'g0': 0, u'cc': [1, 1], u'eqn': [], u'signH': u'[0,2,4,6]', u'con': u'[3,5,5]', u'label': u'2.6-2.0.3-6-6', u'r': 3, u'genus': 2, u'full_label': u'2.24-8.0.2-4-6', u'signature': u'[0,3,6,6]', u'group': u'[6,2]', u'_id': ObjectId('579ba55675d28e4971d3ed7b'), u'group_order': 6})
    cap.insert({u'ndim': 0, u'dim': 0, u'full_auto': u'[10,2]', u'passport_label': u'2.5-1.0.5-5-5.2', u'gen_vectors': [[2, 3, 4, 5, 1], [3, 4, 5, 1, 2], [3, 4, 5, 1, 2]], u'jacobian_decomp': [[2, 1, 2]], u'total_label': u'2.5-1.0.5-5-5.2.1', u'g0': 0, u'cc': [2, 1], u'eqn': [], u'signH': u'[0,2,5,10]', u'con': u'[2,3,3]', u'label': u'2.5-1.0.5-5-5', u'r': 3, u'genus': 2, u'full_label': u'2.10-2.0.2-5-10', u'signature': u'[0,5,5,5]', u'group': u'[5,1]', u'_id': ObjectId('579ba55475d28e4971d3ed78'), u'group_order': 5})
    cap.insert({u'ndim': 0, u'dim': 0, u'full_auto': u'[10,2]', u'passport_label': u'2.5-1.0.5-5-5.3', u'gen_vectors': [[3, 4, 5, 1, 2], [5, 1, 2, 3, 4], [5, 1, 2, 3, 4]], u'jacobian_decomp': [[2, 1, 2]], u'total_label': u'2.5-1.0.5-5-5.3.1', u'g0': 0, u'cc': [3, 1], u'eqn': [], u'signH': u'[0,2,5,10]', u'con': u'[3,5,5]', u'label': u'2.5-1.0.5-5-5', u'r': 3, u'genus': 2, u'full_label': u'2.10-2.0.2-5-10', u'signature': u'[0,5,5,5]', u'group': u'[5,1]', u'_id': ObjectId('579ba55575d28e4971d3ed79'), u'group_order': 5})
    cap.insert({u'ndim': 1, u'dim': 1, u'full_auto': u'[8,3]', u'passport_label': u'2.4-1.0.2-2-4-4.1', u'gen_vectors': [[2, 1, 4, 3], [2, 1, 4, 3], [3, 4, 2, 1], [4, 3, 1, 2]], u'jacobian_decomp': [[2, 1, 2]], u'total_label': u'2.4-1.0.2-2-4-4.1.1', u'g0': 0, u'cc': [1, 1], u'eqn': [], u'signH': u'[0,2,2,2,4]', u'con': u'[2,2,3,4]', u'label': u'2.4-1.0.2-2-4-4', u'r': 4, u'genus': 2, u'full_label': u'2.8-3.0.2-2-2-4', u'signature': u'[0,2,2,4,4]', u'group': u'[4,1]', u'_id': ObjectId('579ba55375d28e4971d3ed76'), u'group_order': 4})
    cap.insert({u'ndim': 0, u'dim': 0, u'full_auto': u'[10,2]', u'passport_label': u'2.5-1.0.5-5-5.1', u'gen_vectors': [[2, 3, 4, 5, 1], [2, 3, 4, 5, 1], [4, 5, 1, 2, 3]], u'jacobian_decomp': [[2, 1, 2]], u'total_label': u'2.5-1.0.5-5-5.1.1', u'g0': 0, u'cc': [1, 1], u'eqn': [], u'signH': u'[0,2,5,10]', u'con': u'[2,2,4]', u'label': u'2.5-1.0.5-5-5', u'r': 3, u'genus': 2, u'full_label': u'2.10-2.0.2-5-10', u'signature': u'[0,5,5,5]', u'group': u'[5,1]', u'_id': ObjectId('579ba55475d28e4971d3ed77'), u'group_order': 5})
