import { useState, useEffect } from 'react';
import Select from 'react-select';

export default function Home() {
  const [abaAtual, setAbaAtual] = useState('gene');
  const [escalaRegiao, setEscalaRegiao] = useState('micro');
  const [opcoes, setOpcoes] = useState([]);
  const [selecionada, setSelecionada] = useState(null);
  const [carregandoLista, setCarregandoLista] = useState(false);
  const [carregandoImagem, setCarregandoImagem] = useState(false);
  const [inputValue, setInputValue] = useState('');

  const [textoListaGenes, setTextoListaGenes] = useState(''); 
  const [imagemCustomizada, setImagemCustomizada] = useState(null);

  const [historicoGene, setHistoricoGene] = useState([]);
  const [historicoOntology, setHistoricoOntology] = useState([]);
  const [historicoGenelist, setHistoricoGenelist] = useState([]);
  const [genesInvalidos, setGenesInvalidos] = useState([]);
  const [mostrarModalGenes, setMostrarModalGenes] = useState(false);

  const [genesDaVia, setGenesDaVia] = useState([]);
  const [mostrarGenesDaVia, setMostrarGenesDaVia] = useState(false);
  const [carregandoGenesVia, setCarregandoGenesVia] = useState(false);

  // Buscando lista de Genes ou Ontologias no R (com Retry automático)
  useEffect(() => {
    setSelecionada(null);
    if (abaAtual === 'genelist') return;

    let montado = true;
    setCarregandoLista(true);
    setOpcoes([]);
    setImagemCustomizada(null);

    setSelecionada(null);

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
          if (formatadas.length > 0) {
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

  // Atualiza os históricos limitando a 3 itens
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

    // 1. Limpa o aviso da pesquisa anterior para não acumular
    setGenesInvalidos([]);
    setMostrarModalGenes(false);

    try {
      const params = new URLSearchParams();
      params.append('gene_string', textoListaGenes);
      params.append('escala', escalaRegiao);

      // 2. CHECAGEM RÁPIDA: Quem existe e quem não existe?
      const resValidacao = await fetch('http://127.0.0.1:33857/validate_genelist', {
        method: 'POST', body: params
      });
      const dadosValidacao = await resValidacao.json();

      // Formata a resposta do R com segurança
      const listaInvalidos = Array.isArray(dadosValidacao.invalidos) ? dadosValidacao.invalidos : dadosValidacao.invalidos?.[0] || [];
      const listaValidos = Array.isArray(dadosValidacao.validos) ? dadosValidacao.validos : dadosValidacao.validos?.[0] || [];

      // Se houver algum erro, salva no Estado (Isso faz a caixa amarela aparecer!)
      if (listaInvalidos.length > 0) {
        setGenesInvalidos(listaInvalidos); 
      }

      // Se TUDO estiver errado, para por aqui e avisa
      if (listaValidos.length === 0) {
        alert("Nenhum dos genes informados foi encontrado na base de dados.");
        setCarregandoImagem(false);
        return; 
      }

      // 3. GERAÇÃO DA IMAGEM
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

  // Historico
  const carregarDoHistorico = (termo) => {
    if (abaAtual !== 'genelist') {
      setSelecionada({ value: termo, label: termo });
      atualizarHistorico(termo, abaAtual);
      setCarregandoImagem(true);
    } else {
      setTextoListaGenes(termo);
    }
  };

  let urlImagem = '';
  if (abaAtual === 'genelist' && imagemCustomizada) {
    urlImagem = imagemCustomizada;
  } else if (selecionada && abaAtual !== 'genelist') {
    urlImagem = abaAtual === 'gene' 
      ? `http://127.0.0.1:33857/plot_brain?gene=${selecionada.value}&escala=${escalaRegiao}`
      : `http://127.0.0.1:33857/plot_ontology?geneset=${selecionada.value}&escala=${escalaRegiao}`;
  }


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
        url = `http://127.0.0.1:33857/data_brain?gene=${selecionada.value}&escala=${escalaRegiao}`;
      } else if (abaAtual === 'ontology') {
        url = `http://127.0.0.1:33857/data_ontology?geneset=${selecionada.value}&escala=${escalaRegiao}`;
      } else if (abaAtual === 'genelist') {
        url = 'http://127.0.0.1:33857/data_genelist';
        const params = new URLSearchParams();
        params.append('gene_string', textoListaGenes);
        params.append('escala', escalaRegiao);
        opcoesFetch = { method: 'POST', body: params };
      }

      const response = await fetch(url, opcoesFetch);
      const blob = await response.blob();
      const objUrl = URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = objUrl;
      a.download = `Dados_${abaAtual}_${escalaRegiao}.csv`;
      a.click();
      URL.revokeObjectURL(objUrl);
    } catch (err) {
      alert("Erro ao baixar os dados da tabela.");
    }
  };

  // Histórico ativo atual
  const historicoAtivo = abaAtual === 'gene' ? historicoGene : abaAtual === 'ontology' ? historicoOntology : historicoGenelist;

  const buscarGenesDaVia = async (nomeOntologia) => {
    if (!nomeOntologia) return;
    
    setCarregandoGenesVia(true);
    setMostrarGenesDaVia(true);
    setGenesDaVia([]); 
    
    try {
      // encodeURIComponent protege nomes complexos na hora de virar URL
      const urlFetch = `http://127.0.0.1:33857/genes_da_via?geneset=${encodeURIComponent(nomeOntologia)}`;
      
      const res = await fetch(urlFetch);
      const dados = await res.json();

      let listaLimpa = [];
      if (dados.genes && Array.isArray(dados.genes)) {
        listaLimpa = dados.genes;
      } else if (Array.isArray(dados)) {
        listaLimpa = dados; 
      } else if (dados.genes && typeof dados.genes === 'string') {
        listaLimpa = [dados.genes];
      }

      setGenesDaVia(listaLimpa);

    } catch (err) {
      console.error("Erro ao buscar genes da via:", err);
      setGenesDaVia(["Erro ao carregar os genes."]);
    } finally {
      setCarregandoGenesVia(false);
    }
  };

  return (
    <section className="max-w-7xl mx-auto px-12 py-12">
      <div className="grid grid-cols-1 lg:grid-cols-12 gap-12">
        
        {/* DIVISÃO ENTRE ABAS MICRO - MACRO REGIÃO */}
        <div className="lg:col-span-4 bg-[#f2f4f4] rounded-xl p-8 border border-zinc-100 h-fit">

          <div className="flex bg-zinc-300 p-1 rounded-lg mb-4">
            {['micro', 'macro'].map((escala) => (
              <button 
                key={escala}
                onClick={() => {
                  setEscalaRegiao(escala);
                  setCarregandoImagem(true);
                }}
                className={`flex-1 py-2 text-[11px] font-black uppercase tracking-widest rounded-md transition-all ${escalaRegiao === escala ? 'bg-[#58614c] text-white shadow-md' : 'text-zinc-500 hover:text-zinc-700'}`}
              >
                {escala === 'micro' ? 'MicroRegions' : 'MacroRegions'}
              </button>
            ))}
          </div>

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
            
            {/* CABEÇALHO COM O BOTÃO DE EXEMPLO (Lado a lado) */}
            <div className="flex justify-between items-center">
              <label className="block text-[10px] font-bold uppercase tracking-widest text-[#58614c]">
                {abaAtual === 'gene' ? 'Search Identifier' : abaAtual === 'ontology' ? 'Search Ontology Term' : 'Input List'}
              </label>

              {/* SÓ MOSTRA O BOTÃO SE ESTIVER NA ABA GENELIST */}
              {abaAtual === 'genelist' && (
                <button
                  onClick={() => setTextoListaGenes("SOX10, TSPAN6, SCYL3, GABRB3, GAD1, GAD2, SLC32A1")}
                  className="text-[10px] text-zinc-500 hover:text-[#58614c] font-bold underline transition-colors"
                >
                  Use example data
                </button>
              )}
            </div>

            {abaAtual !== 'genelist' ? (
              <div className="flex flex-col gap-3">
                <Select
                  value={selecionada}
                  onChange={(opt) => { 
                    setSelecionada(opt); 
                    atualizarHistorico(opt.value, abaAtual);
                    setCarregandoImagem(true); 
                    setMostrarGenesDaVia(false); // Esconde a lista ao trocar de via
                  }}
                  onInputChange={(val) => setInputValue(val)}
                  options={opcoes.filter(o => o.label.toLowerCase().includes(inputValue.toLowerCase())).slice(0, 100)}
                  filterOption={null}
                  styles={{ menu: (base) => ({ ...base, zIndex: 9999 }) }}
                  isLoading={carregandoLista}
                  placeholder="Digite para buscar..."
                  className="text-sm"
                />
                
                {/* BOTÃO E CAIXA DE GENES DA VIA (Só aparece na aba Ontology e se algo estiver selecionado) */}
                {abaAtual === 'ontology' && selecionada && (
                  <div className="flex flex-col gap-2">
                    <button
                      onClick={() => mostrarGenesDaVia ? setMostrarGenesDaVia(false) : buscarGenesDaVia(selecionada.value)}
                      className="text-[10px] font-bold uppercase tracking-widest text-[#58614c] bg-[#e4e9ea] hover:bg-zinc-300 py-2 rounded-lg transition-colors w-full"
                    >
                      {mostrarGenesDaVia ? 'Esconder Genes da Via' : '👁️ Ver Genes desta Via'}
                    </button>

                    {mostrarGenesDaVia && (
                      <div className="bg-white p-3 rounded-lg border border-zinc-200 shadow-inner animate-fade-in">
                        {carregandoGenesVia ? (
                          <p className="text-xs text-zinc-400 text-center py-4">Carregando genes...</p>
                        ) : (
                          <div className="flex flex-col gap-2">
                            <div className="flex justify-between items-center">
                              <span className="text-[10px] font-bold text-zinc-500">{genesDaVia.length} genes encontrados no mapa:</span>
                              <button
                                onClick={() => {
                                  navigator.clipboard.writeText(genesDaVia.join('\n'));
                                  alert('Genes copiados!');
                                }}
                                className="text-[10px] bg-zinc-100 hover:bg-zinc-200 text-zinc-700 px-2 py-1 rounded"
                              >
                                Copiar
                              </button>
                            </div>
                            <textarea
                              readOnly
                              value={genesDaVia.join('\n')}
                              className="w-full bg-zinc-50 text-zinc-600 font-mono text-[10px] p-2 rounded border border-zinc-200 outline-none resize-y h-24"
                            />
                          </div>
                        )}
                      </div>
                    )}
                  </div>
                )}
              </div>
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
                      title={termo}
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
          
            {/* 1. COLE O TRECHO DO AVISO EXATAMENTE AQUI (NO TOPO DA COLUNA) */}
        {genesInvalidos.length > 0 && (
          <div className="bg-[#fff8e6] border border-[#f0dca5] rounded-xl p-5 flex flex-col gap-3 shadow-sm">
            <div className="flex justify-between items-center">
              <p className="text-sm text-[#8c6d1f] font-medium">
                ⚠️ <strong>{genesInvalidos.length} gene(s)</strong> não foram encontrados no banco e foram ignorados.
              </p>
              <div className="flex gap-2">
                <button
                  onClick={() => {
                    navigator.clipboard.writeText(genesInvalidos.join('\n'));
                    alert('Genes copiados para a área de transferência!');
                  }}
                  className="text-[11px] font-bold uppercase tracking-widest bg-white border border-[#f0dca5] text-[#8c6d1f] px-4 py-2 rounded-lg hover:bg-[#fcf4dc] transition-colors"
                >
                  Copiar Genes
                </button>
          
          <button
            onClick={() => setMostrarModalGenes(!mostrarModalGenes)}
            className="text-[11px] font-bold uppercase tracking-widest bg-white border border-[#f0dca5] text-[#8c6d1f] px-4 py-2 rounded-lg hover:bg-[#fcf4dc] transition-colors"
          >
            {mostrarModalGenes ? 'Esconder Lista' : 'Ver Lista'}
          </button>
        </div>
      </div>

      {mostrarModalGenes && (
        <textarea
          readOnly
          value={genesInvalidos.join('\n')}
          className="w-full mt-2 bg-white text-[#8c6d1f] font-mono text-xs p-3 rounded-lg border border-[#f0dca5] outline-none resize-y h-32 cursor-text focus:border-[#8c6d1f]"
        />
      )}
    </div>
  )}

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