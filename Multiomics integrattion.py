
import sys
import networkx as nx
from networkx import Graph as NXGraph
import matplotlib.pyplot as plt
import statistics
import collections
import io
import re
import pandas as pd
import numpy as np
import pickle
import math
import csv
import os
import pickle
from rdflib.extras.external_graph_libs import *
from rdflib import Graph, Literal, URIRef, Namespace, XSD, RDF, RDFS
from rdflib.extras.external_graph_libs import rdflib_to_networkx_graph
from rdflib.namespace import RDF, RDFS
from SPARQLWrapper import SPARQLWrapper, JSON, CSV
from sklearn.preprocessing import StandardScaler

prefixes = '''
PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>
PREFIX wikibase: <http://wikiba.se/ontology#>
PREFIX wdt: <http://www.wikidata.org/prop/direct/>
PREFIX wd: <http://www.wikidata.org/entity/>
PREFIX vg: <http://biohackathon.org/resource/vg#>
PREFIX uniprotkb: <http://purl.uniprot.org/uniprot/>
PREFIX uberon: <http://purl.obolibrary.org/obo/uo#>
PREFIX taxon: <http://purl.uniprot.org/taxonomy/>
PREFIX sp: <http://spinrdf.org/sp#>
PREFIX skos: <http://www.w3.org/2004/02/skos/core#>
PREFIX sio: <http://semanticscience.org/resource/>
PREFIX sh: <http://www.w3.org/ns/shacl#>
PREFIX schema: <http://schema.org/>
PREFIX sachem: <http://bioinfo.uochb.cas.cz/rdf/v1.0/sachem#>
PREFIX rh: <http://rdf.rhea-db.org/>
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX pubmed: <http://rdf.ncbi.nlm.nih.gov/pubmed/>
PREFIX ps: <http://www.wikidata.org/prop/statement/>
PREFIX pq: <http://www.wikidata.org/prop/qualifier/>
PREFIX patent: <http://data.epo.org/linked-data/def/patent/>
PREFIX p: <http://www.wikidata.org/prop/>
PREFIX owl: <http://www.w3.org/2002/07/owl#>
PREFIX orthodbGroup: <http://purl.orthodb.org/odbgroup/>
PREFIX orthodb: <http://purl.orthodb.org/>
PREFIX orth: <http://purl.org/net/orth#>
PREFIX obo: <http://purl.obolibrary.org/obo/>
PREFIX np: <http://nextprot.org/rdf#>
PREFIX nextprot: <http://nextprot.org/rdf/entry/>
PREFIX mnx: <https://rdf.metanetx.org/schema/>
PREFIX mnet: <https://rdf.metanetx.org/mnet/>
PREFIX mesh: <http://id.nlm.nih.gov/mesh/>
PREFIX lscr: <http://purl.org/lscr#>
PREFIX lipidmaps: <https://www.lipidmaps.org/rdf/>
PREFIX keywords: <http://purl.uniprot.org/keywords/>
PREFIX insdcschema: <http://ddbj.nig.ac.jp/ontologies/nucleotide/>
PREFIX insdc: <http://identifiers.org/insdc/>
PREFIX identifiers: <http://identifiers.org/>
PREFIX glyconnect: <https://purl.org/glyconnect/>
PREFIX glycan: <http://purl.jp/bio/12/glyco/glycan#>
PREFIX genex: <http://purl.org/genex#>
PREFIX foaf: <http://xmlns.com/foaf/0.1/>
PREFIX eunisSpecies: <http://eunis.eea.europa.eu/rdf/species-schema.rdf#>
PREFIX ensembltranscript: <http://rdf.ebi.ac.uk/resource/ensembl.transcript/>
PREFIX ensemblterms: <http://rdf.ebi.ac.uk/terms/ensembl/>
PREFIX ensemblprotein: <http://rdf.ebi.ac.uk/resource/ensembl.protein/>
PREFIX ensemblexon: <http://rdf.ebi.ac.uk/resource/ensembl.exon/>
PREFIX ensembl: <http://rdf.ebi.ac.uk/resource/ensembl/>
PREFIX ec: <http://purl.uniprot.org/enzyme/>
PREFIX dcterms: <http://purl.org/dc/terms/>
PREFIX dc: <http://purl.org/dc/terms/>
PREFIX chebislash: <http://purl.obolibrary.org/obo/chebi/>
PREFIX chebihash: <http://purl.obolibrary.org/obo/chebi#>
PREFIX cco: <http://rdf.ebi.ac.uk/terms/chembl#>
PREFIX busco: <http://busco.ezlab.org/schema#>
PREFIX bibo: <http://purl.org/ontology/bibo/>
PREFIX allie: <http://allie.dbcls.jp/>
PREFIX SWISSLIPID: <https://swisslipids.org/rdf/SLM_>
PREFIX GO: <http://purl.obolibrary.org/obo/GO_>
PREFIX ECO: <http://purl.obolibrary.org/obo/ECO_>
PREFIX CHEBI: <http://purl.obolibrary.org/obo/CHEBI_>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX up: <http://purl.uniprot.org/core/>
PREFIX faldo: <http://biohackathon.org/resource/faldo#>
'''

query = '''
SELECT ?protein ?transcript ?ensprotein ?gene
WHERE
{
  ?protein rdfs:seeAlso ?transcript .
  ?protein a up:Protein .
  ?protein up:reviewed true .
  ?protein up:organism taxon:9606 .
  ?transcript a up:Transcript_Resource .
  ?transcript up:translatedTo ?ensprotein .
  ?transcript up:transcribedFrom ?gene .
}
'''

q=f'''
{prefixes}
{query}
'''
class KnowledgeGraph():
    def __init__(self, sparql_endpoint='http://sparql.uniprot.org/sparql'):
        self.sparql_endpoint = sparql_endpoint
    def query_sparql(self, q, output_format='json'):
        sparql = SPARQLWrapper(self.sparql_endpoint)
        sparql.setQuery(q)
        sparql.setReturnFormat(output_format)
        results = sparql.query()
        res = results.convert()
        return res
    def query_to_pandas(self, q):
        res = self.query_sparql(q, output_format='csv')
        df = pd.read_csv(io.BytesIO(res), sep=",")
        return df
    def query_to_graph(self, q):
        res = self.query_sparql(q, output_format='json')  # Change the output_format to 'json'
        # RDF graph setup
        rdf_graph = Graph()
        ens_namespace = Namespace("http://rdf.ebi.ac.uk/resource/ensembl/")
        transcript_namespace = Namespace("http://rdf.ebi.ac.uk/resource/ensembl.transcript/")
        protein_namespace = Namespace("http://rdf.ebi.ac.uk/resource/ensembl.protein/")
        uniprot_namespace = Namespace("http://purl.uniprot.org/uniprot/")

        rdf_graph.bind('ens', ens_namespace)
        rdf_graph.bind('transcript', transcript_namespace)
        rdf_graph.bind('protein', protein_namespace)
        rdf_graph.bind('uniprot', uniprot_namespace)

        # Process bindings
        for binding in res["results"]["bindings"]:
            transcript_uri = URIRef(binding["transcript"]["value"])
            ensprotein_uri = URIRef(binding["ensprotein"]["value"])
            gene_uri = URIRef(binding["gene"]["value"].rsplit('.', 1)[0])
            protein_uri = URIRef(binding["protein"]["value"])

            rdf_graph.add((transcript_uri, RDF.type, transcript_namespace.Transcript))
            rdf_graph.add((ensprotein_uri, RDF.type, protein_namespace.Protein))
            rdf_graph.add((gene_uri, RDF.type, ens_namespace.Gene))
            rdf_graph.add((protein_uri, RDF.type, uniprot_namespace.Protein))

            rdf_graph.add((transcript_uri, transcript_namespace.translateTo, ensprotein_uri))
            rdf_graph.add((gene_uri, ens_namespace.transcribeTo, transcript_uri))
            rdf_graph.add((ensprotein_uri, protein_namespace.isEquivalentTo, protein_uri))
        # Print or serialize the RDF graph
        rdf_graph.serialize(destination='/content/drive/MyDrive/KAUST Master/Research/Heatstroke/r.ttl', format='turtle')
        return rdf_graph

def Knowledge_mirna(r,miRNA_DB, knowledge_graph):
    miRNA_DB['gene'] = miRNA_DB['gene'].apply(add_prefix_en)
    r['gene'] = r['gene'].apply(gene_mod)
    merged_df_mirna_trans = pd.merge(r, miRNA_DB, on='gene', how='left')
    merged_df_mirna_trans['miRNA'] = merged_df_mirna_trans['miRNA'].apply(add_prefix_miRNA)

    mirna_namespace = Namespace("http://rdf.ebi.ac.uk/resource/ensembl.miRNA/")
    knowledge_graph.bind('mirna', mirna_namespace)
    for index, row in merged_df_mirna_trans.iterrows():
        if pd.notna(row['miRNA']) and  pd.notna(row['transcript']):
            subject = URIRef(row['miRNA'])
            obj = URIRef(row['transcript'])
            knowledge_graph.add((subject, RDF.type, mirna_namespace.Gene))
            knowledge_graph.add((subject, mirna_namespace.regulates, obj))
    return knowledge_graph

def Knowledge_mirna_csv(r,miRNA_DB):
    r['gene'] = r['gene'].apply(remove_urls)
    r['gene'] = r['gene'].apply(gene_mod)
    merged_df_mirna_trans = pd.merge(r, miRNA_DB, on='gene', how='left')
    merged_df_mirna_trans['gene'] = merged_df_mirna_trans['gene'].apply(add_prefix_en)
    merged_df_mirna_trans['miRNA'] = merged_df_mirna_trans['miRNA'].apply(add_prefix_miRNA)
    merged_df_mirna_trans.rename(columns={'miRNA':'mirna'}, inplace=True)
    return merged_df_mirna_trans

def add_prefix_en(value):
    if pd.notna(value):
        return 'http://rdf.ebi.ac.uk/resource/ensembl/' + str(value)
    else:
        return value

def add_prefix_t(value):
    if pd.notna(value):  # Check if value is not NaN
        return 'http://rdf.ebi.ac.uk/resource/ensembl.transcript/' + str(value)
    else:
        return value

def add_prefix_uni(value):
    if pd.notna(value):  # Check if value is not NaN
        return 'http://purl.uniprot.org/uniprot/' + str(value)
    else:
        return value

def add_prefix_miRNA(value):
    if pd.notna(value):  # Check if value is not NaN
        return 'http://rdf.ebi.ac.uk/resource/ensembl.miRNA/' + str(value)
    else:
        return value

def gene_mod(value):
       return value.rsplit('.', 1)[0]

def get_uuid(df, file_name, col_name):
    # Filter DataFrame based on file_name
    result_df = df[df['file_name'] == file_name]
    # Check if any matching rows are found
    if not result_df.empty:
        # Extract and return the UUID
        return result_df[col_name].iloc[0]
    else:
        return None

def convert_miRNA_list(mirna_list):
    new_list = [miRNA.split('-5p')[0].split('-3p')[0] if '-5p' in miRNA else
                (miRNA.split('-3p')[0] if '-3p' in miRNA else miRNA) for miRNA in mirna_list]
    return new_list

def check_if_node_exists(graph, node_uri):
    node_ref = URIRef(node_uri)
    return (node_ref, None, None) in graph

def map_proteins(prot_db, prot_table):
    prot = pd.merge(prot_table, prot_db, on='AGID', how='left')
    prot = prot[['protein','protein_expression']]
    prot = prot.dropna(subset=['protein_expression'])
    prot = prot.drop_duplicates(subset=['protein'])
    prot.rename(columns={"gene_id": "gene"}, inplace=True)
    return prot

def process_strings(input_list):
    processed_list = []
    for string in input_list:
        if pd.isna(string):
            processed_list.append(np.nan)
        else:
            parts = string.split('-')
            result = '-'.join(parts[:3])
            processed_list.append(result)
    return processed_list

def iterate_and_read_files(folder_path, omic, disease, target_extension=None): ###1
    omic2=omic+'_expression'
    directory_path = os.path.join(folder_path, omic2, disease)
    if not os.path.isdir(directory_path):
        print(f"The path '{directory_path}' is not a directory.")
        return
    dataframes_dict = {}
    for filename in os.listdir(directory_path):
        file_path = os.path.join(directory_path, filename)
        if os.path.isfile(file_path):
            if target_extension and filename.endswith(target_extension):
                if omic2 == 'gene_expression':
                    file_content = pd.read_csv(file_path, sep='\t', skiprows=[0,2,3,4,5])
                else:
                    file_content = pd.read_csv(file_path, sep='\t')
                dataframes_dict[filename] = file_content   
    return dataframes_dict

def create_rdf_graph(dataframe, omic, UUID, prot_db):
    rdf_graph = Graph()
    ns="http://rdf.ebi.ac.uk/resource/sample/"
    sample_uri = Namespace(ns)
    sample_uri2 = URIRef(ns + UUID)
    rdf_graph.add((sample_uri2, RDF.type, URIRef(omic)))
    rdf_graph.bind('sample', sample_uri)
    rdf_graph.bind('xsd', XSD)
    if omic == 'protein_expression':
        dataframe= map_proteins(prot_db, dataframe)
        uniprot_namespace = Namespace("http://purl.uniprot.org/uniprot/")
        rdf_graph.bind('uniprot', uniprot_namespace)
        #dataframe['gene_id']=URIRef(add_prefix_en(dataframe['gene_id']))
    if omic == 'gene_expression':
        dataframe['gene_id'] = dataframe['gene_id'].str.split('/')
        ens_namespace = Namespace("http://rdf.ebi.ac.uk/resource/ensembl/")       
        rdf_graph.bind('ens', ens_namespace) 
    if omic == 'mirna_expression':
        mirna_namespace = Namespace("http://rdf.ebi.ac.uk/resource/ensembl.miRNA/")
        rdf_graph.bind('mirna', mirna_namespace)
    #rdf_graph.add((sample_uri, RDF.type, RDFS.Resource))
    #rdf_graph.add((custom_ns, custom_ns.hasSample, sample_uri))
    for index, row in dataframe.iterrows():
        if omic == 'gene_expression':
              if pd.notna(row['gene_id']) and  pd.notna(row['tpm_unstranded']):    
                subject = URIRef(add_prefix_en(row['gene_id'][0]))
                #print(subject)
                rdf_graph.add((subject, RDF.type, ens_namespace.Gene))
                rdf_graph.add((sample_uri2, ens_namespace.hasGene, subject))
                obj = Literal(row['tpm_unstranded'], datatype=XSD.double)
                predicate = sample_uri.hasExpression
                rdf_graph.add((subject, predicate, obj))
        elif omic == 'protein_expression':
              if pd.notna(row['protein']) and pd.notna(row[omic]):
                subject = URIRef(add_prefix_uni(row['protein']))
                rdf_graph.add((subject, RDF.type, uniprot_namespace.Protein))
                rdf_graph.add((sample_uri2, uniprot_namespace.hasProtein, subject))
                #rdf_graph.add((sample_uri2, URIRef(sample_uri + 'hasProtein'), subject))
                obj = Literal(row[omic], datatype=XSD.double)
                predicate = sample_uri.hasExpression
                rdf_graph.add((subject, predicate, obj))
        elif omic == 'mirna_expression':
              if pd.notna(row['miRNA_ID']) and pd.notna(row['reads_per_million_miRNA_mapped']):
                subject = URIRef(add_prefix_miRNA(row['miRNA_ID']))
                rdf_graph.add((subject, RDF.type, mirna_namespace.miRNA))
                rdf_graph.add((sample_uri2, mirna_namespace.hasmiRNA, subject))
                #rdf_graph.add((sample_uri2, URIRef(sample_uri + 'hasmiRNA'), subject))
                obj = Literal(row['reads_per_million_miRNA_mapped'], datatype=XSD.double)
                predicate = sample_uri.hasExpression
                rdf_graph.add((subject, predicate, obj))               
    return rdf_graph

def create_rdf_graph_dict(dataframes_dict, df, prot_db):  ####2
    #df = pd.read_csv('/encrypted/e3008/Yang_tcga/metadata/metadata.tsv', sep='\t')
    graphs_dict = {}
    # Iterate through each item in the dictionary
    for filename, dataframe in dataframes_dict.items():
        UUID = get_uuid(df, filename, 'submitter_id_x')
        omic = get_uuid(df, filename, 'type')       
        rdf_graph = create_rdf_graph(dataframe, omic, UUID, prot_db)
        #rdf_graph.serialize(destination=f'/home/castilmg/multiomics/output_{omic}_{UUID}.ttl', format='turtle')
        graphs_dict[filename] = rdf_graph
    return graphs_dict

def create_rdf_graph_study(graphs_dict, df, merged_df, disease):  
    ns = "http://rdf.ebi.ac.uk/resource/study/"
    custom_ns = Namespace(ns)
    study_uri = URIRef(ns + disease)
    rdf_graph = Graph()
    rdf_graph.bind('study', ns)
    rdf_graph.add((study_uri, RDF.type, custom_ns.Study))
    
    for filename, graph in graphs_dict.items():
        UUID = get_uuid(df, filename, 'submitter_id')
        barcode = '-'.join(UUID.split('-')[:3]) 
        omic = get_uuid(df, filename, 'type') 
        ns_sample = Namespace("http://rdf.ebi.ac.uk/resource/sample/")
        rdf_graph.bind('sample', ns_sample)
        sample_uri = URIRef(ns_sample + UUID)
        rdf_graph.add((sample_uri, RDF.type, URIRef(omic)))
        rdf_graph.add((study_uri, custom_ns.hasSample, sample_uri))
        index_of_value = merged_df.index.get_loc(merged_df[merged_df['bcr_patient_barcode'] == barcode].index[0])
        row = merged_df.iloc[index_of_value]
        gender = row['gender']
        age = row['age_at_index']
        project = row['project']
        #type_c = row['primary_pathology_histological_type']
        stage = row['ajcc_pathologic_stage']
        #rdf_graph.add((sample_uri, custom_ns.isOmic, Literal(omic)))
        if pd.notna(gender):
            #rdf_graph.add((sample_uri, URIRef(custom_ns + 'gender'), Literal(gender)))
            rdf_graph.add((sample_uri, ns_sample.isGender, Literal(gender)))
        if pd.notna(age):
            #rdf_graph.add((sample_uri, URIRef(custom_ns + 'age'), Literal(age)))
            rdf_graph.add((sample_uri, ns_sample.hasAge, Literal(age)))
        if pd.notna(project):
            #rdf_graph.add((sample_uri, URIRef(custom_ns + 'projet'), Literal(project)))
            rdf_graph.add((sample_uri, ns_sample.isProject, Literal(project)))
        #if pd.notna(type_c):
            #rdf_graph.add((sample_uri, URIRef(custom_ns + 'type'), Literal(type_c)))
        #    rdf_graph.add((sample_uri, ns_sample.isType, Literal(type_c)))
        if pd.notna(stage):
            #rdf_graph.add((sample_uri, URIRef(custom_ns + 'stage'), Literal(stage)))
            rdf_graph.add((sample_uri, ns_sample.isStage, Literal(stage)))
    
    return rdf_graph

def create_sentence(row,omic):
    if omic=='protein':
        if pd.notna(row['gene']) and pd.notna(row['transcript']):
            return f"The sample {row['Sample']} has protein expression of {row['protein']} with expression level {row['protein_Expression']}, translated from transcript {row['transcript']}, transcribed from gene {row['gene']}."
        else:
            return f"The sample {row['Sample']} has expression level of protein {row['protein']} with expression level {row['protein_Expression']}." 
    elif omic == 'gene':
        if pd.notna(row['protein']) and pd.notna(row['transcript']) and pd.notna(row['mirna']):
            return f"The sample {row['Sample']} has gene expression of {row['gene']} with expression level {row['gene_Expression']}, transcribed to transcript {row['transcript']}, regulated by mirna {row['mirna']}, translated to protein {row['protein']}."
        else:
            return f"The sample {row['Sample']} has expression level of gene {row['gene']} with expression level {row['gene_Expression']}." 
    elif omic == 'mirna':
        if pd.notna(row['gene']) and pd.notna(row['transcript']):
            return f"The sample {row['Sample']} has mirna expression of {row['mirna']} with expression level {row['mirna_Expression']}, that regulates the transcript {row['transcript']} that translates to protein {row['protein']} and is trasncribed from {row['gene']}."
        else:
            return f"The sample {row['Sample']} has expression level of mirna {row['mirna']} with expression level {row['mirna_Expression']}." 
    elif omic == 'multi2':
        if pd.notna(row['gene_Expression']) and pd.notna(row['transcript']) and pd.notna(row['protein_Expression']) and pd.notna(row['mirna_Expression']):
            return f"The sample {row['Sample']} has gene expression of {row['gene']} with expression level {row['gene_Expression']}, transcribed to transcript {row['transcript']}, regulated by mirna {row['mirna']} with expression level {row['mirna_Expression']}, and translated to protein {row['protein']} with expression level {row['protein_Expression']}."
    elif omic == 'multi':
        sentence = " "
        if pd.notna(row['gene']) and pd.notna(row['gene_Expression']):
            s = f"The sample {row['Sample']} has gene expression of {row['gene']} with expression level {row['gene_Expression']}"
            sentence = sentence + s
            if pd.notna(row['transcript']):
                s = f" transcribed to transcript {row['transcript']}"
                sentence = sentence + s 
            if pd.notna(row['protein']):
                s = f" translated to protein {row['protein']}" 
                sentence = sentence + s
                if pd.notna(row['protein_Expression']):
                    s = f" with expression level {row['protein_Expression']}"
                    sentence = sentence + s
            if pd.notna(row['mirna']): 
                s = f" regulated by mirna {row['mirna']}"
                sentence = sentence + s
                if pd.notna(row['mirna_Expression']):
                    s = f" with expression level {row['mirna_Expression']}"
                    sentence = sentence + s
        elif pd.notna(row['protein']) and pd.notna(row['protein_Expression']):
            s = f"The sample {row['Sample']} has protein expression of {row['protein']} with expression level {row['protein_Expression']}"
            sentence = sentence + s
            if pd.notna(row['transcript']):
                s = f" translated from transcript {row['transcript']}"
                sentence = sentence + s 
            if pd.notna(row['gene']):
                s = f" transcribed from gene {row['gene']}"
                sentence = sentence + s
                if pd.notna(row['gene_Expression']):
                    s = f" with expression level {row['gene_Expression']}"
                    sentence = sentence + s
            if pd.notna(row['mirna']): 
                s = f" regulated by mirna {row['mirna']}"
                sentence = sentence + s
                if pd.notna(row['mirna_Expression']):
                    s = f" with expression level {row['mirna_Expression']}"
                    sentence = sentence + s
        elif pd.notna(row['mirna']) and pd.notna(row['mirna_Expression']):
            s = f"The sample {row['Sample']} has mirna expression of {row['mirna']} with expression level {row['mirna_Expression']}"
            sentence = sentence + s
            if pd.notna(row['transcript']):
                s = f" that regulated the transcript {row['transcript']}"
                sentence = sentence + s 
            if pd.notna(row['gene']):
                s = f" transcribed from gene {row['gene']}"
                sentence = sentence + s
                if pd.notna(row['gene_Expression']):
                    s = f" with expression level {row['gene_Expression']}"
                    sentence = sentence + s
            if pd.notna(row['protein']):
                s = f" translated to protein {row['protein']}" 
                sentence = sentence + s
                if pd.notna(row['protein_Expression']):
                    s = f" with expression level {row['protein_Expression']}"
                    sentence = sentence + s
        s = f"the patient is {row['Gender']} with age {row['Age']} the cancer is {row['Type']}"
        sentence = sentence + s
        return sentence

            

def remove_urls(cell):
    if isinstance(cell, str):
        urls_to_remove = [
            'http://rdf.ebi.ac.uk/resource/sample/',
            'http://rdf.ebi.ac.uk/resource/study/',
            'http://rdf.ebi.ac.uk/resource/ensembl.miRNA/',
            'http://rdf.ebi.ac.uk/resource/ensembl/',
            'http://purl.uniprot.org/uniprot/',
            'http://rdf.ebi.ac.uk/resource/ensembl.transcript/',
            'http://rdf.ebi.ac.uk/resource/ensembl.protein/'
        ]
        for url in urls_to_remove:
            cell = cell.replace(url, '')
    return cell



def query_to_df(sample_graph, study_graph, knowledge_graph, omic):
    if omic == 'protein':
        sparql_query_sample = """
        PREFIX sample: <http://rdf.ebi.ac.uk/resource/sample/>
        PREFIX uniprot: <http://purl.uniprot.org/uniprot/>

        SELECT ?sample ?protein ?expression
        WHERE {
            ?sample a <protein_expression> ;
                    uniprot:hasProtein ?protein .

            ?protein a uniprot:Protein ;
                    sample:hasExpression ?expression .
        }
        """
    elif omic == 'gene':
        sparql_query_sample = """
        PREFIX sample: <http://rdf.ebi.ac.uk/resource/sample/>
        PREFIX ens: <http://rdf.ebi.ac.uk/resource/ensembl/>

        SELECT ?sample ?gene ?expression
        WHERE {
            ?sample a <gene_expression> ;
                    ens:hasGene ?gene .

            ?gene a ens:Gene ;
                    sample:hasExpression ?expression .
        }
        """
    elif omic == 'mirna':
        sparql_query_sample = """
        PREFIX sample: <http://rdf.ebi.ac.uk/resource/sample/>
        PREFIX mirna: <http://rdf.ebi.ac.uk/resource/ensembl.miRNA/>

        SELECT ?sample ?mirna ?expression
        WHERE {
            ?sample a <mirna_expression> ;
                    mirna:hasmiRNA ?mirna .

            ?mirna a mirna:miRNA ;
                    sample:hasExpression ?expression .
        }
        """

    sparql_query_study = """
            PREFIX sample: <http://rdf.ebi.ac.uk/resource/sample/>
            PREFIX study: <http://rdf.ebi.ac.uk/resource/study/>
            SELECT ?study ?sample ?age ?gender ?type
            WHERE {
            ?study a study:Study .
            ?study study:hasSample ?sample .
            ?sample sample:hasAge ?age .
            ?sample sample:isGender ?gender .
            ?sample sample:isType ?type .
            }
            """
    sparql_query_kg = """PREFIX ens: <http://rdf.ebi.ac.uk/resource/ensembl/>
            PREFIX protein: <http://rdf.ebi.ac.uk/resource/ensembl.protein/>
            PREFIX transcript: <http://rdf.ebi.ac.uk/resource/ensembl.transcript/>
            PREFIX uniprot: <http://purl.uniprot.org/uniprot/>
            PREFIX mirna: <http://rdf.ebi.ac.uk/resource/ensembl.miRNA/>

            SELECT ?gene ?transcript ?ensprotein ?uniprot ?mirna
            WHERE {
              ?mirna  mirna:regulates ?transcript .
              ?gene ens:hasTranscript ?transcript .
              ?transcript ens:hasProtein ?ensprotein .
              ?ensprotein protein:isEquivalentTo ?uniprot.
            }"""
    # Execute SPARQL query
    sample = sample_graph.query(sparql_query_sample)
    columns = ["Sample", omic, omic+"_Expression"]
    data = [(result["sample"], result[omic], result["expression"]) for result in sample]
    df_sample = pd.DataFrame(data, columns=columns)
    # Execute SPARQL query
    study = study_graph.query(sparql_query_study)
    columns = ["Study","Sample", "Age", "Gender", "Type"]
    data = [(result["study"], result["sample"], result["age"], result["gender"], result["type"]) for result in study]
    print(data)
    df_study = pd.DataFrame(data, columns=columns)
    # Execute SPARQL query
    kg = knowledge_graph.query(sparql_query_kg)
    columns = ["gene","transcript", "ens_protein", "protein", "mirna"]
    data = [(result["gene"], result["transcript"], result["ensprotein"], result["uniprot"], result["mirna"]) for result in kg]
    df_kg = pd.DataFrame(data, columns=columns)

    merged_df_sample_kg = pd.merge(df_sample, df_kg, on=omic, how='left')
    merged_study = pd.merge(merged_df_sample_kg, df_study, on='Sample', how='left')
    merged_df_sample_kg = merged_df_sample_kg.map(remove_urls)
    merged_study = merged_study.map(remove_urls)
    return merged_df_sample_kg, merged_study


def query_to_df(sample_graph, study_graph, omic):
    if omic == 'protein':
        sparql_query_sample = """
        PREFIX sample: <http://rdf.ebi.ac.uk/resource/sample/>
        PREFIX uniprot: <http://purl.uniprot.org/uniprot/>

        SELECT ?sample ?protein ?expression
        WHERE {
            ?sample a <protein_expression> ;
                    uniprot:hasProtein ?protein .

            ?protein a uniprot:Protein ;
                    sample:hasExpression ?expression .
        }
        """
    elif omic == 'gene':
        sparql_query_sample = """
        PREFIX sample: <http://rdf.ebi.ac.uk/resource/sample/>
        PREFIX ens: <http://rdf.ebi.ac.uk/resource/ensembl/>

        SELECT ?sample ?gene ?expression
        WHERE {
            ?sample a <gene_expression> ;
                    ens:hasGene ?gene .

            ?gene a ens:Gene ;
                    sample:hasExpression ?expression .
        }
        """
    elif omic == 'mirna':
        sparql_query_sample = """
        PREFIX sample: <http://rdf.ebi.ac.uk/resource/sample/>
        PREFIX mirna: <http://rdf.ebi.ac.uk/resource/ensembl.miRNA/>

        SELECT ?sample ?mirna ?expression
        WHERE {
            ?sample a <mirna_expression> ;
                    mirna:hasmiRNA ?mirna .

            ?mirna a mirna:miRNA ;
                    sample:hasExpression ?expression .
        }
        """

    sparql_query_study = """
            PREFIX sample: <http://rdf.ebi.ac.uk/resource/sample/>
            PREFIX study: <http://rdf.ebi.ac.uk/resource/study/>
            SELECT ?study ?sample ?age ?gender ?type
            WHERE {
            ?study a study:Study .
            ?study study:hasSample ?sample .
            ?sample sample:hasAge ?age .
            ?sample sample:isGender ?gender .
            ?sample sample:isType ?type .
            }
            """
    sparql_query_kg = """PREFIX ens: <http://rdf.ebi.ac.uk/resource/ensembl/>
            PREFIX protein: <http://rdf.ebi.ac.uk/resource/ensembl.protein/>
            PREFIX transcript: <http://rdf.ebi.ac.uk/resource/ensembl.transcript/>
            PREFIX uniprot: <http://purl.uniprot.org/uniprot/>
            PREFIX mirna: <http://rdf.ebi.ac.uk/resource/ensembl.miRNA/>

            SELECT ?gene ?transcript ?ensprotein ?uniprot ?mirna
            WHERE {
              ?mirna  mirna:regulates ?transcript .
              ?gene ens:hasTranscript ?transcript .
              ?transcript ens:hasProtein ?ensprotein .
              ?ensprotein protein:isEquivalentTo ?uniprot.
            }"""
    
    
    # Execute SPARQL query
    sample = sample_graph.query(sparql_query_sample)
    columns = ["Sample", omic, omic+"_Expression"]
    data = [(result["sample"], result[omic], result["expression"]) for result in sample]
    df_sample = pd.DataFrame(data, columns=columns)
    # Execute SPARQL query
    study = study_graph.query(sparql_query_study)
    columns = ["Study","Sample", "Age", "Gender", "Type"]
    data = [(result["study"], result["sample"], result["age"], result["gender"], result["type"]) for result in study]
    df_study = pd.DataFrame(data, columns=columns)   
    
    merged_study_sample = pd.merge(df_sample, df_study, on='Sample', how='left') 
    merged_study_sample = merged_study_sample.map(remove_urls)
    
    if omic == 'gene':
        merged_study_sample[omic]=merged_study_sample[omic].apply(gene_mod)

    return merged_study_sample

def split_dict(d, n):
    """
    Split a dictionary into n smaller dictionaries.

    Parameters:
    - d (dict): The dictionary to be split.
    - n (int): The number of smaller dictionaries to split the original dictionary into.

    Returns:
    - list: A list of smaller dictionaries.
    """
    if n <= 0:
        raise ValueError("Number of splits (n) must be greater than 0.")
    if len(d) < n:
        raise ValueError("Number of splits (n) cannot exceed the length of the dictionary.")

    # Calculate the approximate size of each smaller dictionary
    split_size = len(d) // n

    # Split the dictionary into smaller dictionaries
    smaller_dicts = []
    start_idx = 0
    for i in range(n):
        end_idx = start_idx + split_size
        smaller_dict = dict(list(d.items())[start_idx:end_idx])
        smaller_dicts.append(smaller_dict)
        start_idx = end_idx

    # Distribute the remaining items to the last smaller dictionary
    if start_idx < len(d):
        last_smaller_dict = dict(list(d.items())[start_idx:])
        smaller_dicts[-1].update(last_smaller_dict)

    return smaller_dicts

def create_study_df(graph, omic):
    if omic == 'protein':
        sparql_query_sample = """
        PREFIX sample: <http://rdf.ebi.ac.uk/resource/sample/>
        PREFIX uniprot: <http://purl.uniprot.org/uniprot/>

        SELECT ?sample ?protein ?expression
        WHERE {
            ?sample a <protein_expression> ;
                    uniprot:hasProtein ?protein .

            ?protein a uniprot:Protein ;
                    sample:hasExpression ?expression .
        }
        """
    elif omic == 'gene':
        sparql_query_sample = """
        PREFIX sample: <http://rdf.ebi.ac.uk/resource/sample/>
        PREFIX ens: <http://rdf.ebi.ac.uk/resource/ensembl/>

        SELECT ?sample ?gene ?expression
        WHERE {
            ?sample a <gene_expression> ;
                    ens:hasGene ?gene .

            ?gene a ens:Gene ;
                    sample:hasExpression ?expression .
        }
        """
    elif omic == 'mirna':
        sparql_query_sample = """
        PREFIX sample: <http://rdf.ebi.ac.uk/resource/sample/>
        PREFIX mirna: <http://rdf.ebi.ac.uk/resource/ensembl.miRNA/>

        SELECT ?sample ?mirna ?expression
        WHERE {
            ?sample a <mirna_expression> ;
                    mirna:hasmiRNA ?mirna .

            ?mirna a mirna:miRNA ;
                    sample:hasExpression ?expression .
        }
        """
        
    sam_study=pd.DataFrame()
    for filename, sample_graph in graph.items():
        sample = sample_graph.query(sparql_query_sample)
        columns = ["Sample", omic, omic+"_Expression"]
        data = [(result["sample"], result[omic], result["expression"]) for result in sample]
        df_sample = pd.DataFrame(data, columns=columns)
        sam_study = pd.concat([sam_study, df_sample], ignore_index=True)
    sam_study = sam_study.map(remove_urls)
    return sam_study
    

disease = 'TCGA-BRCA'
prot_DB = pd.read_csv('/home/castilmg/multiomics/TCGA_antibodies_descriptions.gencode.v36_2.tsv',sep='\t')

patient_df = pd.read_csv(F"/home/castilmg/multiomics/clinical_data_{disease}.csv")
#patient_df = patient_df[['bcr_patient_barcode','age_at_index', 'gender', 'project', 'primary_pathology_histological_type', 'stage_event_pathologic_stage']]
patient_df['merge_column'] = process_strings(patient_df['bcr_patient_barcode'])
#patient_df.to_csv('/home/castilmg/multiomics/patient_df.tsv', sep='\t', index=False, header=True)

meta = pd.read_csv('/encrypted/e3008/Yang_tcga/metadata/metadata.tsv', sep='\t')
meta = meta[meta['cases'] == disease]
meta = meta[(meta['data_format'] != 'BCR XML') & (meta['data_format'] != 'PDF') & (meta['data_format'] != 'BCR Biotab') & (meta['data_format'] != 'BCR OMF XML')]
meta = meta[(meta['type'] == 'mirna_expression') | (meta['type'] == 'gene_expression') | (meta['type'] == 'protein_expression')]
meta['submitter_id'].fillna(meta['platform'], inplace=True)
meta['merge_column'] = process_strings(meta['submitter_id'])
meta = meta[['UUID','id','data_format','cases','file_name','data_category','type','submitter_id','merge_column']]
meta.to_csv('/home/castilmg/multiomics/meta.csv', index=False, header=True)

merged_df = pd.merge(meta, patient_df, on='merge_column', how='left')

#merged_df = pd.read_csv('/home/castilmg/multiomics/merged_meta_df.csv')

merged_df.to_csv(f"/home/castilmg/multiomics/merged_df_{disease}.tsv", sep='\t', index=False, header=True)

folder_path='/encrypted/e3008/Yang_tcga/tcga/'
extension_to_read = '.tsv'  # Specify the extension you want to read, or set to None

# PROTEINS

omic = 'protein'
prot_dic=iterate_and_read_files(folder_path, omic, disease, extension_to_read)
prot_graph=create_rdf_graph_dict(prot_dic, merged_df, prot_DB)
with open(f"/home/castilmg/multiomics/prot_graph_{disease}.pkl", 'wb') as file:
    pickle.dump(prot_graph, file)
print("done proteins")

# GENES

omic = 'gene'
trans=iterate_and_read_files(folder_path, omic, disease, extension_to_read)
smaller_dicts = split_dict(trans, 4)
c=0
for sm in smaller_dicts:
    trans_graph=create_rdf_graph_dict(sm, merged_df, prot_DB)
    c=c+1
    with open(f"/home/castilmg/multiomics/trans_graph_{disease}_{c}.pkl", 'wb') as file:
        pickle.dump(trans_graph, file)
    print(f"done genes_{c}") 

# miRNAs
omic = 'mirna'
extension_to_read = 'mirnas.quantification.txt'
miRNAs=iterate_and_read_files(folder_path, omic, disease, extension_to_read)
miRNA_graph=create_rdf_graph_dict(miRNAs, merged_df, prot_DB)
with open(f"/home/castilmg/multiomics/miRNA_graph_{disease}.pkl", 'wb') as file:
    pickle.dump(miRNA_graph, file)
print("done mirnas")

# Multiomics
with open(f"/home/castilmg/multiomics/prot_graph_{disease}.pkl", 'rb') as file:
    prot_graph = pickle.load(file)
with open(f"/home/castilmg/multiomics/miRNA_graph_{disease}.pkl", 'rb') as file:
    miRNA_graph = pickle.load(file)
multiomics_graph_dic = {**prot_graph, **miRNA_graph}
multiomics_study_graph=create_rdf_graph_study(multiomics_graph_dic, meta, patient_df, disease)

multiomics_study_graph.serialize(destination=f'/home/castilmg/multiomics/output_study_{disease}_1.ttl', format='turtle')
with open(f"/home/castilmg/multiomics/multiomics_graph_dic_{disease}_1.pkl", 'wb') as file:
    pickle.dump(multiomics_graph_dic, file)    
print("done Multiomics 1")


with open(f"/home/castilmg/multiomics/trans_graph_{disease}_1.pkl", 'rb') as file:
    multiomics_graph_dic = pickle.load(file)
multiomics_study_graph=create_rdf_graph_study(multiomics_graph_dic, meta, patient_df, disease)
multiomics_study_graph.serialize(destination=f'/home/castilmg/multiomics/output_study_{disease}_2.ttl', format='turtle')
print("done Multiomics 2")

with open(f"/home/castilmg/multiomics/trans_graph_{disease}_2.pkl", 'rb') as file:
    multiomics_graph_dic = pickle.load(file)   
multiomics_study_graph=create_rdf_graph_study(multiomics_graph_dic, meta, patient_df, disease)
multiomics_study_graph.serialize(destination=f'/home/castilmg/multiomics/output_study_{disease}_3.ttl', format='turtle')
print("done Multiomics 3")
    
with open(f"/home/castilmg/multiomics/trans_graph_{disease}_3.pkl", 'rb') as file:
    multiomics_graph_dic = pickle.load(file)
multiomics_study_graph=create_rdf_graph_study(multiomics_graph_dic, meta, patient_df, disease)
multiomics_study_graph.serialize(destination=f'/home/castilmg/multiomics/output_study_{disease}_4.ttl', format='turtle')
print("done Multiomics 4")

multiomics_study_graph1 = Graph()
multiomics_study_graph1.parse(f"/home/castilmg/multiomics/output_study_{disease}_1.ttl", format='turtle')
multiomics_study_graph2 = Graph()
multiomics_study_graph2.parse(f"/home/castilmg/multiomics/output_study_{disease}_2.ttl", format='turtle')
multiomics_study_graph3 = Graph()
multiomics_study_graph3.parse(f"/home/castilmg/multiomics/output_study_{disease}_3.ttl", format='turtle')
multiomics_study_graph4 = Graph()
multiomics_study_graph4.parse(f"/home/castilmg/multiomics/output_study_{disease}_4.ttl", format='turtle')

multiomics_study_graph = multiomics_study_graph1 + multiomics_study_graph2 + multiomics_study_graph3 + multiomics_study_graph4

multiomics_study_graph.serialize(destination=f'/home/castilmg/multiomics/output_study_{disease}.ttl', format='turtle')
print("done Multiomics")
    


#multiomics_graph_dic = {**prot_graph, **trans_graph1, **miRNA_graph}
#print("done Multiomics")

# Save
#with open(f"/home/castilmg/multiomics/multiomics_graph_dic_{disease}_1.pkl", 'wb') as file:
#    pickle.dump(multiomics_graph_dic, file)

#multiomics_graph_dic = {**trans_graph2, **trans_graph3}
#print("done Multiomics 2")

# Save
#with open(f"/home/castilmg/multiomics/multiomics_graph_dic_{disease}_2.pkl", 'wb') as file:
#    pickle.dump(multiomics_graph_dic, file)    
    
# Study
with open(f"/home/castilmg/multiomics/multiomics_graph_dic_{disease}.pkl", 'rb') as file:
    multiomics_graph_dic = pickle.load(file)

multiomics_study_graph=create_rdf_graph_study(multiomics_graph_dic, meta, patient_df, disease)
print("DONE Graph")
multiomics_study_graph.serialize(destination=f'/home/castilmg/multiomics/output_study_{disease}_1.ttl', format='turtle')
print("DONE")

knowledge_graph = Graph()
knowledge_graph.parse('/home/castilmg/multiomics/Kownlege_graph.ttl', format='turtle')

miRNA_DB = pd.read_csv('/home/castilmg/multiomics/miRTarBase.csv')
miRNA_DB.loc[:, 'miRNA'] = convert_miRNA_list(miRNA_DB['miRNA'])
miRNA_DB.rename(columns={'ENSEMBL':'gene','Target Gene':'SYMBOL'}, inplace=True)
miRNA_DB = miRNA_DB[['miRNA','gene']]

kg = KnowledgeGraph()
df_kg = kg.query_to_pandas(q)
df_kg.columns=['protein','transcript','ensprotein', 'gene']
print(df_kg.head())

#knowledge_graph = Knowledge_mirna(df_kg,miRNA_DB, knowledge_graph)
#knowledge_graph.serialize(destination='/home/castilmg/multiomics/Kownlege_graph.ttl', format='turtle')

multiomics_study_graph = Graph()
multiomics_study_graph.parse(f"/home/castilmg/multiomics/output_study_{disease}.ttl", format='turtle')

#knowledge_graph = Graph()
#knowledge_graph.parse('/home/castilmg/multiomics/Kownlege_graph.ttl', format='turtle')

df_kg2 = Knowledge_mirna_csv(df_kg,miRNA_DB)
df_kg2=df_kg2.map(remove_urls)
print(df_kg2.head())

# PROTEIN
with open(f"/home/castilmg/multiomics/prot_graph_{disease}.pkl", 'rb') as file:
    prot_graph = pickle.load(file)
omic = 'protein'
protein_sam_study = create_study_df(prot_graph, omic)

# GENES
with open(f"/home/castilmg/multiomics/trans_graph_{disease}.pkl", 'rb') as file:
    trans_graph = pickle.load(file)
omic = 'gene'
gene_sam_study = create_study_df(trans_graph, omic)
gene_sam_study['gene'] = gene_sam_study['gene'].apply(gene_mod)
# miRNA
with open(f"/home/castilmg/multiomics/miRNA_graph_{disease}.pkl", 'rb') as file:
    miRNA_graph = pickle.load(file)
omic = 'mirna'
mirna_sam_study = create_study_df(miRNA_graph, omic)

# Merge multiomics

prot = pd.merge(protein_sam_study, df_kg, on='protein', how='outer')
#print(prot)
gene_ = pd.merge(gene_sam_study, prot, on=['Sample','gene'], how='inner')
#print(gene_)
multiomics_df = pd.merge(mirna_sam_study, gene_, on=['mirna','Sample'], how='outer')
#print(multiomics_df)


# Merge Study

sparql_query_study = """
    PREFIX sample: <http://rdf.ebi.ac.uk/resource/sample/>
    PREFIX study: <http://rdf.ebi.ac.uk/resource/study/>
    SELECT ?study ?sample ?age ?gender ?type
    WHERE {
    ?study a study:Study .
    ?study study:hasSample ?sample .
    ?sample sample:hasAge ?age .
    ?sample sample:isGender ?gender .
    ?sample sample:isType ?type .
    }
    """


study = multiomics_study_graph.query(sparql_query_study)
columns = ["Study","Sample", "Age", "Gender", "Type"]
data = [(result["study"], result["sample"], result["age"], result["gender"], result["type"]) for result in study]
df_study = pd.DataFrame(data, columns=columns)
df_study = df_study.map(remove_urls)

merged_multiomics_study = pd.merge(multiomics_df, df_study, on='Sample', how='left')
merged_multiomics_study['sentence']=merged_multiomics_study.apply(create_sentence, axis=1, args=(omic,))
columns_to_standardize = ['mirna_Expression', 'gene_Expression', 'protein_Expression']


#NORMALIZE THE OMICS
merged_multiomics_study['mirna_Expression']=(merged_multiomics_study['mirna_Expression']-merged_multiomics_study['mirna_Expression'].min())/(merged_multiomics_study['mirna_Expression'].max()-merged_multiomics_study['mirna_Expression'].min())#normalize

merged_multiomics_study['gene_Expression']=(merged_multiomics_study['gene_Expression']-merged_multiomics_study['gene_Expression'].min())/(merged_multiomics_study['gene_Expression'].max()-merged_multiomics_study['gene_Expression'].min())#normalize

merged_multiomics_study['protein_Expression']=(merged_multiomics_study['protein_Expression']-merged_multiomics_study['protein_Expression'].min())/(merged_multiomics_study['protein_Expression'].max()-merged_multiomics_study['protein_Expression'].min())#normalize
#scaler = StandardScaler()# Use StandardScaler
#merged_multiomics_study[columns_to_standardize] = scaler.fit_transform(merged_multiomics_study[columns_to_standardize])
omic='multi'
merged_multiomics_study['sentence'] = merged_multiomics_study.apply(create_sentence, axis=1, args=(omic,))
merged_multiomics_study.to_csv(f"/home/castilmg/multiomics/multiomics_study_df_{disease}.csv")


#MERGE WITH STUDY
# PROTEIN
with open('/home/castilmg/multiomics/prot_graph.pkl', 'rb') as file:
    prot_graph = pickle.load(file)
protein_df_kg_study=pd.DataFrame()
for filename, p_graph in prot_graph.items():
    omic = 'protein'
    merged_df_sample_kg, merged_study = query_to_df(p_graph, multiomics_study_graph, knowledge_graph, omic)
    protein_df_kg_study = pd.concat([protein_df_kg_study, merged_study], ignore_index=True)
protein_df_kg_study['sentence'] = protein_df_kg_study.apply(create_sentence, axis=1, args=(omic,))
protein_df_kg_study.to_csv('/home/castilmg/multiomics/protein_df_kg_study.csv', index=False, header=True)
file = f"/home/castilmg/multiomics/{disease}_{omic}_samples.txt"
with open(file, 'a') as f:
    df_string = protein_df_kg_study['sentence'].to_string(header=False, index=False)
    f.write(df_string)
# GENES
with open('/home/castilmg/multiomics/trans_graph.pkl', 'rb') as file:
    trans_graph = pickle.load(file)
gene_df_kg_study=pd.DataFrame()
for filename, gene_graph in trans_graph.items():
    omic = 'gene'
    merged_df_sample_kg, merged_study = query_to_df(gene_graph, multiomics_study_graph, knowledge_graph, omic)
    gene_df_kg_study = pd.concat([gene_df_kg_study, merged_study], ignore_index=True)
gene_df_kg_study['sentence'] = gene_df_kg_study.apply(create_sentence, axis=1, args=(omic,))
gene_df_kg_study.to_csv('/home/castilmg/multiomics/gene_df_kg_study.csv', index=False, header=True)
file = f"/home/castilmg/multiomics/{disease}_{omic}_samples.txt"
with open(file, 'a') as f:
    df_string = gene_df_kg_study['sentence'].to_string(header=False, index=False)
    f.write(df_string)
# miRNA
with open('/home/castilmg/multiomics/miRNA_graph.pkl', 'rb') as file:
    miRNA_graph = pickle.load(file)
mirna_df_kg_study=pd.DataFrame()
for filename, mi_graph in miRNA_graph.items():
    omic = 'mirna'
    merged_df_sample_kg, merged_study = query_to_df(mi_graph, multiomics_study_graph, knowledge_graph, omic)
    mirna_df_kg_study = pd.concat([mirna_df_kg_study, merged_study], ignore_index=True)
mirna_df_kg_study['sentence'] = mirna_df_kg_study.apply(create_sentence, axis=1, args=(omic,))
mirna_df_kg_study.to_csv('/home/castilmg/multiomics/mirna_df_kg_study.csv', index=False, header=True)
file = f"/home/castilmg/multiomics/{disease}_{omic}_samples.txt"
with open(file, 'a') as f:
    df_string = mirna_df_kg_study['sentence'].to_string(header=False, index=False)
    f.write(df_string)



"""import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader
from sklearn.model_selection import train_test_split
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.preprocessing import MinMaxScaler

# Sample data (replace this with your data)
data = [
    ("The sample has expression level of gene X with expression level 10.", 10),
    ("Another sample has expression level of gene Y with expression level 15.", 15),
    # Add more data samples
]
Data = pd.read_csv('/home/castilmg/multiomics/multiomics_df_kg_study.csv',sep=',')
Data = Data.applymap(remove_urls)
Data_ = Data['sentence'].to_list()
expression_levels = [float(match.group(1)) for sentence in Data_ for match in re.finditer(r'expression level ([\d.-]+)\b', sentence)]

data = [(sentence, expression_level) for sentence, expression_level in zip(Data_, expression_levels)]
# Split data into train and test sets
X = [sample[0] for sample in data]
Y = [sample[1] for sample in data]
X_train, X_test, y_train, y_test = train_test_split(X,Y, test_size=0.2, random_state=42)

# Tokenize and vectorize the text data
vectorizer = CountVectorizer()
X_train_text = vectorizer.fit_transform(X_train)
X_test_text = vectorizer.transform(X_test)

# Normalize the expression levels
scaler = MinMaxScaler()
y_train = scaler.fit_transform(np.array(y_train).reshape(-1, 1))
y_test = scaler.transform(np.array((y_test)).reshape(-1, 1))

# Convert to PyTorch tensors
X_train_text = torch.tensor(X_train_text.toarray(), dtype=torch.float32)
X_test_text = torch.tensor(X_test_text.toarray(), dtype=torch.float32)
y_train = torch.tensor(y_train, dtype=torch.float32)
y_test = torch.tensor(y_test, dtype=torch.float32)

# Define a simple neural network for self-supervised learning
class ExpressionPredictor(nn.Module):
    def __init__(self, input_size, hidden_size):
        super(ExpressionPredictor, self).__init__()
        self.fc = nn.Linear(input_size, hidden_size)
        self.relu = nn.ReLU()
        self.output_layer = nn.Linear(hidden_size, 1)

    def forward(self, x):
        x = self.fc(x)
        x = self.relu(x)
        x = self.output_layer(x)
        return 

# Train the model
input_size = X_train_text.shape[1]
hidden_size = 64  # Adjust as needed
model = ExpressionPredictor(input_size, hidden_size)
criterion = nn.MSELoss()
optimizer = torch.optim.Adam(model.parameters(), lr=0.001)

# Training loop
num_epochs = 10  # Adjust as needed
for epoch in range(num_epochs):
    optimizer.zero_grad()
    outputs = model(X_train_text)
    loss = criterion(outputs, y_train)
    loss.backward()
    optimizer.step()
    print(f"Epoch [{epoch+1}/{num_epochs}], Loss: {loss.item()}")

# Evaluate the model on the test set
with torch.no_grad():
    test_outputs = model(X_test_text)
    test_loss = criterion(test_outputs, y_test)
    print(f"Test Loss: {test_loss.item()}")

"""
"""
import torch
from transformers import GPT2LMHeadModel, GPT2Tokenizer
import random
from torch.utils.data import Dataset
from transformers import Trainer, TrainingArguments
from sklearn.model_selection import train_test_split
from transformers import default_data_collator

data_collator = default_data_collator

def mask_random_words(sentence, mask_symbol='<mask>', mask_probability=0.2):
    words = re.findall(r'\b\w+\b', sentence)  # Extract words from the sentence
    masked_sentence = sentence

    for word in words:
        if random.random() < mask_probability:
            masked_sentence = re.sub(r'\b' + re.escape(word) + r'\b', mask_symbol, masked_sentence)

    return masked_sentence

class TextCompletionDataset(Dataset):
    def __init__(self, input_texts, output_texts, tokenizer, max_length):
        self.input_texts = input_texts
        self.output_texts = output_texts
        self.tokenizer = tokenizer
        self.max_length = max_length

    def __len__(self):
        return len(self.input_texts)

    def __getitem__(self, idx):
        input_text = self.input_texts[idx]
        output_text = self.output_texts[idx]

        # Tokenize input and output texts
        input_ids = self.tokenizer.encode(input_text, max_length=self.max_length, return_tensors="pt", truncation=True, padding='max_length')[0]
        output_ids = self.tokenizer.encode(output_text, max_length=self.max_length, return_tensors="pt", truncation=True, padding='max_length')[0]

        return {"input_ids": input_ids, "labels": output_ids}

# Add a padding token to the tokenizer

Data = pd.read_csv('/home/castilmg/multiomics/multiomics_df_kg_study.csv',sep=',')
Data_ = Data['sentence'].to_list()
Data['mask'] = Data['sentence'].apply(mask_random_words)

X=Data['mask'].to_list()
Y=Data['sentence'].to_list()
x_train_data, x_test_data, y_train_data, y_test_data, = train_test_split(X,Y, test_size=0.2, random_state=42)



model_name = "gpt2"
tokenizer = GPT2Tokenizer.from_pretrained(model_name)
tokenizer.add_special_tokens({'pad_token': '[PAD]'})  # Add a padding token

model = GPT2LMHeadModel.from_pretrained(model_name)
train_dataset = TextCompletionDataset(x_train_data, y_train_data, tokenizer, 256)

training_args = TrainingArguments(
    output_dir="./text_completion_model",
    overwrite_output_dir=True,
    num_train_epochs=1,
    per_device_train_batch_size=64,
)
trainer = Trainer(
    model=model,
    args=training_args,
    data_collator=data_collator,
    train_dataset=train_dataset,
)
model.resize_token_embeddings(len(tokenizer))
assert model.transformer.wte.weight.shape[0] == len(tokenizer)
trainer.train()

input_text_test = "/n".join(x_test_data)
output_text_test = "/n".join(y_test_data)
input_ids = tokenizer.encode(input_text_test, return_tensors="pt")

generated_ids = model.generate(input_ids, max_length=150, num_return_sequences=1, no_repeat_ngram_size=2, top_k=50, top_p=0.95)
generated_text = tokenizer.decode(generated_ids[0], skip_special_tokens=True)

print("Generated Text:" , generated_text)
print("Real text:" , output_text_test)

"""