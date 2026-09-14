export default function About() {
  return (
    <section className="max-w-7xl mx-auto px-6 sm:px-8 lg:px-12 py-16">
      <div className="max-w-4xl mx-auto">
        <p className="text-[10px] font-bold text-zinc-400 uppercase tracking-[0.3em] mb-4">
          About This App
        </p>
        <h1 className="text-6xl font-bold text-[#2d3435] leading-tight tracking-tighter mb-8">
          Human Developmental Brain RNA-Seq Dataset
        </h1>

        <div className="text-xl text-[#596061] leading-relaxed max-w-2xl mb-16 space-y-4">
          <p>
            This application lets you explore gene expression across human
            brain development, using RNA-Seq data from the BrainSpan Atlas.
            It's built for researchers who want to move from a gene, an
            ontology term, or a custom gene list straight to an expression
            pattern, without setting up their own analysis pipeline.
          </p>
          <p>
            The underlying data comes from the{' '}
            <a
              href="https://www.brainspan.org/"
              target="_blank"
              rel="noopener noreferrer"
              className="text-[#2d3435] underline decoration-zinc-300 hover:decoration-[#2d3435] transition-colors"
            >
              BrainSpan Atlas of the Developing Human Brain
            </a>
            , which profiles gene expression across brain regions and
            developmental stages, from early fetal periods through adulthood.
          </p>
        </div>

        <div className="max-w-3xl mb-16">
          <h2 className="text-sm font-bold text-[#2d3435] uppercase tracking-[0.2em] mb-6">
            What you can do here
          </h2>
          <div className="grid grid-cols-2 gap-x-12 border-t border-zinc-200">
            <div className="py-6 border-b border-zinc-200">
              <h3 className="text-lg font-semibold text-[#2d3435] mb-1">
                Search by gene
              </h3>
              <p className="text-base text-[#596061] leading-relaxed">
                Look up a single gene and see how its expression changes
                across brain regions and developmental time points.
              </p>
            </div>
            <div className="py-6 border-b border-l border-zinc-200 pl-12">
              <h3 className="text-lg font-semibold text-[#2d3435] mb-1">
                Score custom gene lists with ssGSEA
              </h3>
              <p className="text-base text-[#596061] leading-relaxed">
                Upload your own list of genes and get single-sample GSEA
                enrichment scores computed against the BrainSpan expression
                data, so you can test whether your set of interest shifts
                across regions or stages.
              </p>
            </div>
            <div className="py-6">
              <h3 className="text-lg font-semibold text-[#2d3435] mb-1">
                Browse by ontology
              </h3>
              <p className="text-base text-[#596061] leading-relaxed">
                Explore expression patterns for genes grouped under
                biological ontology terms, to see how whole processes or
                pathways trend across development.
              </p>
            </div>
            <div className="py-6 border-l border-zinc-200 pl-12">
              <h3 className="text-lg font-semibold text-[#2d3435] mb-1">
                Plot custom log2FC values
              </h3>
              <p className="text-base text-[#596061] leading-relaxed">
                Bring in your own log2 fold-change values from an external
                comparison — a different disease model, treatment, or
                condition — and plot them against the BrainSpan expression
                data. This lets you run analyses along directions the
                dataset wasn't originally built for, not just the
                comparisons already included in it.
              </p>
            </div>
          </div>
        </div>

    {/*}
        <div className="max-w-2xl mb-16">
            <h2 className="text-sm font-bold text-[#2d3435] uppercase tracking-[0.2em] mb-4">
              Citation
            </h2>
            <p className="text-base text-[#596061] leading-relaxed">
              Vamo usar quando tiver paper{' '}
              <a
                href="link do paper publicado"
                target="_blank"
                rel="noopener noreferrer"
                className="text-[#2d3435] underline decoration-zinc-300 hover:decoration-[#2d3435] transition-colors"
              >
                Titulo do Paper
              </a>{' '}
              resto da frase.
            </p>
        </div>
    {*/}

        <div className="max-w-2xl">
          <h2 className="text-sm font-bold text-[#2d3435] uppercase tracking-[0.2em] mb-4">
            Built by
          </h2>
          <p className="text-base text-[#596061] leading-relaxed">
            This application is developed and maintained by the{' '}
            <a
              href="https://mcblab.github.io/"
              target="_blank"
              rel="noopener noreferrer"
              className="text-[#2d3435] underline decoration-zinc-300 hover:decoration-[#2d3435] transition-colors"
            >
              Marques-Coelho Bioinformatics Lab
            </a>
            , it is situated in Instituto do Cérebro - ICE, UFRN, Natal-RN.
          </p>
        </div>
 
      </div>
    </section>
  );
}