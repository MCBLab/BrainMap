import { useState, useEffect } from 'react';
import Select from 'react-select';

export default function Home() {
  const [abaAtual, setAbaAtual] = useState('gene'); 
  const [opcoes, setOpcoes] = useState([]);
  const [selecionada, setSelecionada] = useState(null);
  const [carregandoLista, setCarregandoLista] = useState(false);
  const [carregandoImagem, setCarregandoImagem] = useState(false);
  const [inputValue, setInputValue] = useState('');

  // ESTADOS DA GENE LIST
  const [textoListaGenes, setTextoListaGenes] = useState(''); 
  const [imagemCustomizada, setImagemCustomizada] = useState(null);

  // ESTADOS DE HISTÓRICO (Até 3 itens)
  const [historicoGene, setHistoricoGene] = useState([]);
  const [historicoOntology, setHistoricoOntology] = useState([]);
  const [historicoGenelist, setHistoricoGenelist] = useState([]);

  // Busca lista de Genes ou Ontologias no R (com Retry automático)
  useEffect(() => {
    setSelecionada(null);
    if (abaAtual === 'genelist') return;

    let montado = true;
    setCarregandoLista(true);
    setOpcoes([]);
    setImagemCustomizada(null);

    const endpoint = abaAtual === 'gene' ? '/list_genes' : '/list_ontologies';

    const buscarComRetry = () => {
      fetch(`http://127.0.0.1:33857${endpoint}`)
        .then((res) => {
          if (!res.ok) throw new Error("API ligando...");
          return res.json();
        })
        .then((dados) => {
          if (!montado) return;
          const arrayDados = Array.isArray(dados) ? dados : dados[0] || [];
          const formatadas = arrayDados.map(item => ({ value: item, label: item }));
          
          setOpcoes(formatadas);
          setCarregandoLista(false);
          
          // Se não houver nada selecionado, seleciona o primeiro
          if (!selecionada && formatadas.length > 0) {
            setSelecionada(formatadas[0]);
            atualizarHistorico(formatadas[0].value, abaAtual);
            setCarregandoImagem(true);
          }
        })
        .catch(() => { if (montado) setTimeout(buscarComRetry, 2000); });
    };

    buscarComRetry();
    return () => { montado = false; };
  }, [abaAtual]);

  // FUNÇÃO: Atualiza os históricos limitando a 3 itens sem duplicatas
  const atualizarHistorico = (termo, aba) => {
    if (!termo || termo.trim() === '') return;
    
    if (aba === 'gene') {
      setHistoricoGene(prev => [termo, ...prev.filter(t => t !== termo)].slice(0, 3));
    } else if (aba === 'ontology') {
      setHistoricoOntology(prev => [termo, ...prev.filter(t => t !== termo)].slice(0, 3));
    } else if (aba === 'genelist') {
      setHistoricoGenelist(prev => [termo, ...prev.filter(t => t !== termo)].slice(0, 3));
    }
  };

  // FUNÇÃO: Processar Genelist Nova
  const gerarMapaLista = async () => {
    if (!textoListaGenes.trim()) return; 
    
    setCarregandoImagem(true);
    setSelecionada({ label: "Lista Customizada" }); 
    atualizarHistorico(textoListaGenes, 'genelist');

    try {
      const params = new URLSearchParams();
      params.append('gene_string', textoListaGenes);

      const resposta = await fetch('http://127.0.0.1:33857/plot_genelist', {
        method: 'POST', body: params
      });
      if (!resposta.ok) throw new Error("Erro");

      const blob = await resposta.blob();
      setImagemCustomizada(URL.createObjectURL(blob));
    } catch (erro) {
      console.error(erro);
      alert("Erro ao processar a lista.");
    } finally {
      setCarregandoImagem(false);
    }
  };

  // FUNÇÃO: Carregar do Histórico
  const carregarDoHistorico = (termo) => {
    if (abaAtual !== 'genelist') {
      setSelecionada({ value: termo, label: termo });
      atualizarHistorico(termo, abaAtual);
      setCarregandoImagem(true);
    } else {
      setTextoListaGenes(termo);
    }
  };

  // Define URL da imagem ativa
  let urlImagem = '';
  if (abaAtual === 'genelist' && imagemCustomizada) {
    urlImagem = imagemCustomizada;
  } else if (selecionada && abaAtual !== 'genelist') {
    urlImagem = abaAtual === 'gene' 
      ? `http://127.0.0.1:33857/plot_brain?gene=${selecionada.value}`
      : `http://127.0.0.1:33857/plot_ontology?geneset=${selecionada.value}`;
  }

  // ==========================================
  // FUNÇÕES DE DOWNLOAD (SVG e DADOS)
  // ==========================================
  const baixarSVG = async () => {
    if (!urlImagem) return;
    try {
      const response = await fetch(urlImagem);
      const blob = await response.blob();
      const url = URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url;
      a.download = `Mapa_${abaAtual}.svg`;
      a.click();
      URL.revokeObjectURL(url);
    } catch (err) {
      alert("Erro ao tentar baixar o SVG.");
    }
  };

  const baixarDados = async () => {
    try {
      let url = '';
      let opcoesFetch = {};

      if (abaAtual === 'gene') {
        url = `http://127.0.0.1:33857/data_brain?gene=${selecionada.value}`;
      } else if (abaAtual === 'ontology') {
        url = `http://127.0.0.1:33857/data_ontology?geneset=${selecionada.value}`;
      } else if (abaAtual === 'genelist') {
        url = 'http://127.0.0.1:33857/data_genelist';
        const params = new URLSearchParams();
        params.append('gene_string', textoListaGenes);
        opcoesFetch = { method: 'POST', body: params };
      }

      const response = await fetch(url, opcoesFetch);
      const blob = await response.blob();
      const objUrl = URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = objUrl;
      a.download = `Dados_${abaAtual}.csv`;
      a.click();
      URL.revokeObjectURL(objUrl);
    } catch (err) {
      alert("Erro ao baixar os dados da tabela.");
    }
  };

  // Histórico ativo atual
  const historicoAtivo = abaAtual === 'gene' ? historicoGene : abaAtual === 'ontology' ? historicoOntology : historicoGenelist;

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
                onChange={(opt) => { 
                  setSelecionada(opt); 
                  atualizarHistorico(opt.value, abaAtual);
                  setCarregandoImagem(true); 
                }}
                onInputChange={(val) => setInputValue(val)}
                options={opcoes.filter(o => o.label.toLowerCase().includes(inputValue.toLowerCase())).slice(0, 100)}
                filterOption={null}
                styles={{ menu: (base) => ({ ...base, zIndex: 9999 }) }}
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

          {/* HISTÓRICO RECENTE */}
          {historicoAtivo.length > 0 && (
            <div className="mt-8 pt-6 border-t border-zinc-200">
              <span className="text-[10px] font-bold uppercase tracking-widest text-zinc-400 block mb-3">Recent:</span>
              <div className="flex flex-wrap gap-2">
                {historicoAtivo.map((termo, idx) => (
                  <button 
                    key={idx}
                    onClick={() => carregarDoHistorico(termo)}
                    className="bg-white hover:bg-zinc-100 border border-zinc-200 text-zinc-600 text-[11px] font-medium px-3 py-1.5 rounded-full transition-colors max-w-full truncate"
                    title={termo} // Permite ver a lista completa de genes se passar o mouse por cima
                  >
                    {termo.length > 25 ? termo.substring(0, 25) + '...' : termo}
                  </button>
                ))}
              </div>
            </div>
          )}

        </div>

        {/* COLUNA DIREITA (Visualização e Botões) */}
        <div className="lg:col-span-8 flex flex-col gap-4">
          
          {/* Caixa do Mapa */}
          <div className="bg-white border border-zinc-100 rounded-xl overflow-hidden shadow-sm p-8 flex items-center justify-center min-h-[550px] relative">
            {carregandoImagem && <p className="absolute font-bold text-[#58614c] animate-pulse z-10 bg-white/80 px-4 py-2 rounded">⏳ Processando mapa...</p>}
            {urlImagem && (
              <img 
                src={urlImagem} 
                className="max-w-full h-auto transition-opacity duration-300" 
                style={{ opacity: carregandoImagem ? 0.2 : 1 }}
                onLoad={() => setCarregandoImagem(false)}
                alt="Brain Map"
              />
            )}
          </div>

          {/* Botões de Download */}
          {urlImagem && (
            <div className="flex gap-4 justify-end">
              <button 
                onClick={baixarDados}
                className="flex items-center gap-2 px-6 py-2.5 bg-white border border-zinc-200 text-zinc-700 rounded-lg text-xs font-bold uppercase tracking-widest hover:border-zinc-300 hover:bg-zinc-50 transition-all shadow-sm"
              >
                📊 Download Data (CSV)
              </button>
              <button 
                onClick={baixarSVG}
                className="flex items-center gap-2 px-6 py-2.5 bg-[#58614c] text-white rounded-lg text-xs font-bold uppercase tracking-widest hover:opacity-90 shadow-sm transition-all"
              >
                🖼️ Download Plot (SVG)
              </button>
            </div>
          )}

        </div>

      </div>
    </section>
  );
}