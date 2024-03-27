
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
    
    
#Integrate Knowledge graph with miRNAs
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

#Integrate Knowledge graph with miRNAs return in csv
def Knowledge_mirna_csv(r,miRNA_DB):
    r['gene'] = r['gene'].apply(remove_urls)
    r['gene'] = r['gene'].apply(gene_mod)
    merged_df_mirna_trans = pd.merge(r, miRNA_DB, on='gene', how='left')
    merged_df_mirna_trans['gene'] = merged_df_mirna_trans['gene'].apply(add_prefix_en)
    merged_df_mirna_trans['miRNA'] = merged_df_mirna_trans['miRNA'].apply(add_prefix_miRNA)
    merged_df_mirna_trans.rename(columns={'miRNA':'mirna'}, inplace=True)
    return merged_df_mirna_trans

def convert_miRNA_list(mirna_list):
    new_list = [miRNA.split('-5p')[0].split('-3p')[0] if '-5p' in miRNA else
                (miRNA.split('-3p')[0] if '-3p' in miRNA else miRNA) for miRNA in mirna_list]
    return new_list

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
    split_size = len(d) // n
    smaller_dicts = []
    start_idx = 0
    for i in range(n):
        end_idx = start_idx + split_size
        smaller_dict = dict(list(d.items())[start_idx:end_idx])
        smaller_dicts.append(smaller_dict)
        start_idx = end_idx
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
path= './Multiomics TCGA'
prot_DB = pd.read_csv(path+'/TCGA_antibodies_descriptions.gencode.v36_2.tsv',sep='\t')

merged_df = pd.read_csv(path+'/merged_meta_df.csv')

extension_to_read = '.tsv'  # Specify the extension you want to read, or set to None

# PROTEINS
omic = 'protein'
prot_dic=iterate_and_read_files(path, omic, disease, extension_to_read)
prot_graph=create_rdf_graph_dict(prot_dic, merged_df, prot_DB)
with open(f"{path}/prot_graph_{disease}.pkl", 'wb') as file:
    pickle.dump(prot_graph, file)
print("done proteins")

# GENES
omic = 'gene'
trans=iterate_and_read_files(path, omic, disease, extension_to_read)
trans_graph=create_rdf_graph_dict(trans, merged_df, prot_DB)
with open(f"{path}/trans_graph_{disease}.pkl", 'wb') as file:
    pickle.dump(trans_graph, file)
print(f"done transcripts") 

# miRNAs
omic = 'mirna'
extension_to_read = '.txt'
miRNAs=iterate_and_read_files(path, omic, disease, extension_to_read)
miRNA_graph=create_rdf_graph_dict(miRNAs, merged_df, prot_DB)
with open(f"{path}/miRNA_graph_{disease}.pkl", 'wb') as file:
    pickle.dump(miRNA_graph, file)
print("done mirnas")

# Multiomics
with open(f"{path}/trans_graph_{disease}.pkl", 'rb') as file:
    trans_graph = pickle.load(file)
with open(f"{path}/prot_graph_{disease}.pkl", 'rb') as file:
    prot_graph = pickle.load(file)
with open(f"{path}/miRNA_graph_{disease}.pkl", 'rb') as file:
    miRNA_graph = pickle.load(file)
multiomics_graph_dic = {**trans_graph, **prot_graph, **miRNA_graph}
multiomics_study_graph=create_rdf_graph_study(multiomics_graph_dic, meta, patient_df, disease)
multiomics_study_graph.serialize(destination=f'{path}/output_study_{disease}.ttl', format='turtle')
with open(f"{path}/multiomics_graph_dic_{disease}.pkl", 'wb') as file:
    pickle.dump(multiomics_graph_dic, file)    
print("done Multiomics")

multiomics_study_graph = Graph()
multiomics_study_graph.parse(f"{path}/output_study_{disease}.ttl", format='turtle')

# Study
with open(f"{path}/multiomics_graph_dic_{disease}.pkl", 'rb') as file:
    multiomics_graph_dic = pickle.load(file)

multiomics_study_graph=create_rdf_graph_study(multiomics_graph_dic, meta, patient_df, disease)
multiomics_study_graph.serialize(destination=f'{path}/output_study_{disease}_1.ttl', format='turtle')
print("Done Study Graph")

#Integrate miRNAs to the Knowledge graph
knowledge_graph = Graph()
knowledge_graph.parse(path+'/Kownlege_graph.ttl', format='turtle')

miRNA_DB = pd.read_csv(path+'/miRTarBase.csv')
miRNA_DB.loc[:, 'miRNA'] = convert_miRNA_list(miRNA_DB['miRNA'])
miRNA_DB.rename(columns={'ENSEMBL':'gene','Target Gene':'SYMBOL'}, inplace=True)
miRNA_DB = miRNA_DB[['miRNA','gene']]

kg = KnowledgeGraph()
df_kg = kg.query_to_pandas(q)
df_kg.columns=['protein','transcript','ensprotein', 'gene']

multiomics_study_graph = Graph()
multiomics_study_graph.parse(f"{path}/output_study_{disease}.ttl", format='turtle')

df_kg2 = Knowledge_mirna_csv(df_kg,miRNA_DB)
df_kg2=df_kg2.map(remove_urls)

# Tranform Graph to Dataframes 
# PROTEIN
with open(f"{path}/prot_graph_{disease}.pkl", 'rb') as file:
    prot_graph = pickle.load(file)
omic = 'protein'
protein_sam_study = create_study_df(prot_graph, omic)
# GENES
with open(f"{path}/trans_graph_{disease}.pkl", 'rb') as file:
    trans_graph = pickle.load(file)
omic = 'gene'
gene_sam_study = create_study_df(trans_graph, omic)
gene_sam_study['gene'] = gene_sam_study['gene'].apply(gene_mod)
# miRNA
with open(f"{path}/miRNA_graph_{disease}.pkl", 'rb') as file:
    miRNA_graph = pickle.load(file)
omic = 'mirna'
mirna_sam_study = create_study_df(miRNA_graph, omic)

# Merge individual sample information with Knowledge graph
protein_df = pd.merge(protein_sam_study, df_kg, on='protein', how='outer')
gene_df = pd.merge(gene_sam_study, protein_df, on=['Sample','gene'], how='inner')
multiomics_df = pd.merge(mirna_sam_study, gene_df, on=['mirna','Sample'], how='outer')

# Query Study Graph
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
# Organize the study graph query results in a Dataframe
columns = ["Study","Sample", "Age", "Gender", "Type"]
data = [(result["study"], result["sample"], result["age"], result["gender"], result["type"]) for result in study]
df_study = pd.DataFrame(data, columns=columns)
df_study = df_study.map(remove_urls)

# Merge multiomics dataframe and the Study dataframe
merged_multiomics_study = pd.merge(multiomics_df, df_study, on='Sample', how='left')

#Standardize the Omics
columns_to_standardize = ['mirna_Expression', 'gene_Expression', 'protein_Expression']
merged_multiomics_study['mirna_Expression']=(merged_multiomics_study['mirna_Expression']-merged_multiomics_study['mirna_Expression'].mean())/(merged_multiomics_study['mirna_Expression'].std())#Standardize
merged_multiomics_study['gene_Expression']=(merged_multiomics_study['gene_Expression']-merged_multiomics_study['gene_Expression'].mean())/(merged_multiomics_study['gene_Expression'].std())#Standardize
merged_multiomics_study['protein_Expression']=(merged_multiomics_study['protein_Expression']-merged_multiomics_study['protein_Expression'].mean())/(merged_multiomics_study['protein_Expression'].std())#Standardize

omic='multi'
merged_multiomics_study['sentence'] = merged_multiomics_study.apply(create_sentence, axis=1, args=(omic,))
merged_multiomics_study.to_csv(f"{path}/multiomics_study_df_{disease}.csv")

