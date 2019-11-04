import pandas as pd
import sys
import random
from math import sqrt

class GSEA:
    """
    GSEA performs Gene Set Enrichment Analysis on a given gene set
    """

    def load_data(self, expfile, sampfile, genesets):
        """
        Takes the expression, sample and gene set data, reads them in, and stores them within the GSEA

        Inputs:
            exp_file is the file name for the gene expression data set 
            samp_file is file name for the data set of samples and corresponding labels
            gene_sets is the file name for the list of gene sets
        """
        self.expressions = pd.read_csv(expfile, sep='\t', index_col='SYMBOL')
        self.samples = pd.read_csv(sampfile, sep='\t', header=None, index_col=0)
        self.genesets = {}
        self.genes_ranked = None
        with open(genesets, 'r') as f:
            for line in f.readlines():
                line = line.strip()
                name, url, *genes = line.split('\t')
                genes = filter(lambda g: g in self.expressions.index, genes)
                self.genesets[name] = set(genes)

    def get_expr_means(self, healthy):
        """
        Returns a dataframe containing the mean expression value per gene of patients or control
        
        Inputs:
            patient: bool, whether to return expression values for patients or control
        """
        op = lambda x: x < 0.5 if healthy else x > 0.5
        names = self.samples[op(self.samples[1])].index
        expr = self.expressions[names]
        return expr.mean(axis=1)

    def get_gene_rank_order(self):
        """
        Returns a list of all genes (as strings) in descending logFC order
        """
        healthy_means = self.get_expr_means(True)
        disease_means = self.get_expr_means(False)
        log_fc = disease_means - healthy_means
        log_fc.sort_values(inplace=True, ascending=False)
        return log_fc.index.tolist()

    def get_enrichment_score(self, geneset_name):
        """
        Returns the enrichment score, a float corrected to two decimal places, for a given gene set

        Inputs:
            gene_set is the given gene set
        """
        geneset = self.genesets[geneset_name]
        if self.genes_ranked is None:
            self.genes_ranked = self.get_gene_rank_order()
        es = 0
        cur_sum = 0
        G = len(geneset)
        N = len(self.genes_ranked)

        up_step = sqrt((N-G)/float(G))
        down_step = sqrt(G/float(N-G))

        for gene in self.genes_ranked:
            if gene in geneset:
                cur_sum += up_step
            else:
                cur_sum -= down_step
            if cur_sum > es:
                es = cur_sum

        return es

    def get_p_values(self):
        """
        Returns a dict mapping each geneset name to the p value for whether or not those genes are differentially more expressed
        """
        act_scores = {name: self.get_enrichment_score(name) for name in self.genesets}
        num_greater = {name: 0 for name in self.genesets}
        for i in range(100):
            #permute sample groups
            #random.shuffle(self.genes_ranked)
            for geneset in self.genesets:
                print(geneset)
                score = self.get_enrichment_score(geneset)
                if score > act_scores[geneset]:
                    num_greater[geneset] += 1
        self.genes_ranked = self.get_gene_rank_order()
        return {geneset: num_greater[geneset]/100.0 for geneset in self.genesets}

    def get_sig_sets(self, p):
        """
        Returns the list of significant gene sets (as strings), at a corrected threshold, by name
        Return an empty list if no gene sets are significant

        Inputs:
            p is the threshold epresenting the corrected p-value after Bonferroni
        """
        p_values = self.get_p_values()
        n = len(self.genesets)
        alpha = p/float(n)
        return [geneset for geneset in self.genesets if p_values[geneset] < alpha]
 
def main():

    #check that the file is being properly used
    if (len(sys.argv) !=4):
        print("Please specify an expression file, a sample file, and a kegg file as args.")
        return
        
    #input variables
    expfile = sys.argv[1]
    sampfile = sys.argv[2]
    keggfile = sys.argv[3]

    gsea = GSEA()
    gsea.load_data(expfile, sampfile, keggfile)
    print(gsea.get_sig_sets(0.01))

if __name__ == "__main__":
    main()
