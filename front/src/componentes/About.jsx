export default function About() {
  return (
    <section className="max-w-7xl mx-auto px-12 py-16">
      <div className="max-w-4xl">
        <p className="text-[10px] font-bold text-zinc-400 uppercase tracking-[0.3em] mb-4">
          About This App
        </p>
        <h1 className="text-6xl font-bold text-[#2d3435] leading-tight tracking-tighter mb-8">
          Human Developmental Brain RNA-Seq Dataset.
        </h1>
        <div className="text-xl text-[#596061] leading-relaxed max-w-2xl mb-16 space-y-4">
          <p>
            This application visualizes gene expression data from the Human Developmental Brain RNA-Seq dataset.
          </p>
          <p>
            The data is sourced from the{' '}
            <a 
              href="https://www.brainspan.org/" 
              target="_blank" 
              rel="noopener noreferrer"
              className="text-[#2d3435] underline decoration-zinc-300 hover:decoration-[#2d3435] transition-colors"
            >
              BrainSpan Atlas of the Developing Human Brain
            </a>
            .
          </p>
        </div>
      </div>
    </section>
  );
}