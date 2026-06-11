export default function Header({ paginaAtual, setPaginaAtual }) {
  return (
    <header className="w-full bg-white border-b border-zinc-100 z-50">
      <div className="flex justify-between items-center px-8 py-4 max-w-full mx-auto">
        <div className="text-xl font-bold tracking-tighter text-[#58614c]">VEGABrain</div>
        <nav className="flex gap-x-8">
          <button 
            onClick={() => setPaginaAtual('home')} 
            className={`font-medium text-sm transition-all ${paginaAtual === 'home' ? 'text-[#58614c] border-b-2 border-[#58614c]' : 'text-[#596061] hover:text-[#2d3435]'}`}
          >
            Home
          </button>
          <button 
            onClick={() => setPaginaAtual('about')} 
            className={`font-medium text-sm transition-all ${paginaAtual === 'about' ? 'text-[#58614c] border-b-2 border-[#58614c]' : 'text-[#596061] hover:text-[#2d3435]'}`}
          >
            About
          </button>
        </nav>
      </div>
    </header>
  );
}