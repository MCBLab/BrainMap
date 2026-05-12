import { useState } from 'react';
import Header from './componentes/Header';
import Home from './componentes/Home';
import About from './componentes/About';
import './App.css';

function App() {
  const [paginaAtual, setPaginaAtual] = useState('home');

  return (
    <div className="flex flex-col min-h-screen font-sans bg-white">
      {/* O Header recebe a função para mudar de página */}
      <Header paginaAtual={paginaAtual} setPaginaAtual={setPaginaAtual} />
      
      <main className="flex-grow">
        {/* Renderização Condicional */}
        {paginaAtual === 'home' && <Home />}
        {paginaAtual === 'about' && <About />}
      </main>

      {/* Footer (Opcional, vindo do modelo do Stitch) */}
      <footer className="w-full bg-white border-t border-zinc-100 px-12 py-6">
        <div className="text-[10px] font-medium uppercase tracking-[0.2em] text-[#596061]">
          © 2024 ExpressionBrainMap. Genomic Precision Framework.
        </div>
      </footer>
    </div>
  );
}

export default App;