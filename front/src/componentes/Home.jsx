import { useState, useEffect } from 'react';
import Select from 'react-select';

export default function Home() {
  const [abaAtual, setAbaAtual] = useState('gene'); 
  const [opcoes, setOpcoes] = useState([]);
  const [selecionada, setSelecionada] = useState(null);
  const [carregandoLista, setCarregandoLista] = useState(false);
  const [carregandoImagem, setCarregandoImagem] = useState(false);
  const [inputValue, setInputValue] = useState('');

  // ESTADOS DA NOVA ABA 'GENE LIST'
  const [textoListaGenes, setTextoListaGenes] = useState(''); 
  const [imagemCustomizada, setImagemCustomizada] = useState(null);

  // Busca dados no Plumber R
  useEffect(() => {
    if (abaAtual === 'genelist') return;

    setCarregandoLista(true);
    setOpcoes([]);
    setSelecionada(null);
    setInputValue('');
    setImagemCustomizada(null); // Limpa imagem customizada se trocar de aba

    const endpoint = abaAtual === 'gene' ? '/list_genes' : '/list_ontologies';

    fetch(`http://127.0.0.1:33857${endpoint}`)
      .then((res) => res.json())
      .then((dados) => {
        // Garantindo compatibilidade com o retorno do R
        const arrayDados = Array.isArray(dados) ? dados : dados[0] || [];
        const formatadas = arrayDados.map(item => ({ value: item, label: item }));
        
        setOpcoes(formatadas);
        setCarregandoLista(false);
        
        // A SUA SACADA AQUI: Seleciona o primeiro item automaticamente!
        if (formatadas.length > 0) {
          setSelecionada(formatadas[0]);
          setCarregandoImagem(true);
        }
      })
      .catch(() => setCarregandoLista(false));
  }, [abaAtual]);

  // FUNÇÃO DA GENE LIST: Envia a lista digitada para o R
  const gerarMapaLista = async () => {
    if (!textoListaGenes.trim()) return; 
    
    setCarregandoImagem(true);
    setSelecionada({ label: "Lista Customizada" }); 

    try {
      const params = new URLSearchParams();
      params.append('gene_string', textoListaGenes);

      const resposta = await fetch('http://127.0.0.1:33857/plot_genelist', {
        method: 'POST',
        body: params
      });

      if (!resposta.ok) throw new Error("Erro no servidor");

      // Transforma a resposta em imagem
      const blob = await resposta.blob();
      const urlTemporaria = URL.createObjectURL(blob);
      setImagemCustomizada(urlTemporaria);

    } catch (erro) {
      console.error(erro);
      alert("Erro ao processar a lista. Verifique os genes e tente novamente.");
    } finally {
      setCarregandoImagem(false);
    }
  };

  // Define qual imagem mostrar
  let urlImagem = 'https://raw.githubusercontent.com/google/brave-ui/main/assets/brain.png';
  if (abaAtual === 'genelist' && imagemCustomizada) {
    urlImagem = imagemCustomizada; // Imagem gerada pelo POST
  } else if (selecionada && abaAtual !== 'genelist') {
    // Imagem gerada pelo GET (sua lógica original)
    urlImagem = abaAtual === 'gene' 
      ? `http://127.0.0.1:33857/plot_brain?gene=${selecionada.value}`
      : `http://127.0.0.1:33857/plot_ontology?geneset=${selecionada.value}`;
  }

  return (
    <section className="max-w-7xl mx-auto px-12 py-12">
      <div className="grid grid-cols-1 lg:grid-cols-12 gap-12">
        
        {/* COLUNA ESQUERDA (Análise) */}
        <div className="lg:col-span-4 bg-[#f2f4f4] rounded-xl p-8 border border-zinc-100 h-fit">
          <div className="flex bg-[#e4e9ea] p-1 rounded-lg mb-6">
            {['gene', 'ontology', 'genelist'].map((tab) => (
              <button 
                key={tab}
                onClick={() => setAbaAtual(tab)}
                className={`flex-1 py-2 text-[10px] font-bold uppercase tracking-widest rounded-md transition-all ${abaAtual === tab ? 'bg-white text-[#58614c]' : 'text-[#596061]'}`}
              >
                {tab === 'genelist' ? 'Gene List' : tab}
              </button>
            ))}
          </div>

          <div className="space-y-4">
            <label className="block text-[10px] font-bold uppercase tracking-widest text-[#58614c]">
              {abaAtual === 'gene' ? 'Search Identifier' : abaAtual === 'ontology' ? 'Search Ontology Term' : 'Input List'}
            </label>

            {abaAtual !== 'genelist' ? (
              <Select
                value={selecionada}
                onChange={(opt) => { setSelecionada(opt); setCarregandoImagem(true); }}
                onInputChange={(val) => setInputValue(val)}
                options={opcoes
                  .filter(o => o.label.toLowerCase().includes(inputValue.toLowerCase()))
                  .slice(0, 100)
                }
                filterOption={null}
                styles={{
                  menu: (base) => ({ ...base, zIndex: 9999 }),
                  option: (base, state) => ({ 
                    ...base, color: '#2d3435', backgroundColor: state.isFocused ? '#f2f4f4' : 'white' 
                  })
                }}
                isLoading={carregandoLista}
                placeholder="Digite para buscar..."
                className="text-sm"
              />
            ) : (
              <div>
                <textarea 
                  rows="6" 
                  value={textoListaGenes}
                  onChange={(e) => setTextoListaGenes(e.target.value)}
                  placeholder="Cole sua lista de genes aqui..." 
                  className="w-full px-4 py-3 bg-white border border-transparent focus:border-[#58614c] rounded-lg text-sm outline-none resize-none"
                ></textarea>
                
                <button 
                  onClick={gerarMapaLista}
                  className="w-full mt-4 bg-[#58614c] text-white py-3 rounded-lg font-bold text-xs uppercase tracking-widest hover:opacity-90 active:scale-[0.98] transition-all"
                >
                  Renderizar Mapa
                </button>
              </div>
            )}
          </div>
        </div>

        {/* COLUNA DIREITA (Visualização) */}
        <div className="lg:col-span-8 bg-white border border-zinc-100 rounded-xl overflow-hidden shadow-sm">
          <div className="p-8 flex items-center justify-center min-h-[550px] relative">
            {carregandoImagem && <p className="absolute font-bold text-[#58614c] animate-pulse z-10 bg-white/80 px-4 py-2 rounded">⏳ Processando mapa...</p>}
            
            <img 
              src={urlImagem} 
              className="max-w-full h-auto transition-opacity duration-300" 
              style={{ opacity: carregandoImagem ? 0.2 : 1 }}
              onLoad={() => setCarregandoImagem(false)}
              alt="Brain Map"
            />
          </div>
        </div>

      </div>
    </section>
  );
}