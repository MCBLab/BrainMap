import { useState, useEffect } from 'react';
import Select from 'react-select'; // Importando a nossa nova caixa de busca
import './App.css';

function App() {
  const [modo, setModo] = useState('ontologia'); 
  const [opcoes, setOpcoes] = useState([]);
  
  // Agora o selecionado guarda um objeto {value: '...', label: '...'} por causa do react-select
  const [selecionada, setSelecionada] = useState(null); 
  
  const [carregandoLista, setCarregandoLista] = useState(false);
  
  // NOVO ESTADO: Controla se a imagem está sendo gerada pelo R
  const [carregandoImagem, setCarregandoImagem] = useState(false);

  useEffect(() => {
    setCarregandoLista(true);
    setOpcoes([]);
    setSelecionada(null);

    const endpoint = modo === 'gene' ? '/list_genes' : '/list_ontologies';

    fetch(`http://127.0.0.1:33857${endpoint}`)
      .then((resposta) => resposta.json())
      .then((dados) => {
        // O react-select exige que os dados tenham esse formato de label e value
        const opcoesFormatadas = dados.map(item => ({ value: item, label: item }));
        
        setOpcoes(opcoesFormatadas);
        if (opcoesFormatadas.length > 0) {
          setSelecionada(opcoesFormatadas[0]);
          setCarregandoImagem(true); // Começa a carregar a primeira imagem
        }
        setCarregandoLista(false);
      })
      .catch((erro) => {
        console.error("Erro ao buscar dados da API:", erro);
        setCarregandoLista(false);
      });
  }, [modo]);

  // Monta a URL baseada no 'value' do item selecionado
  const urlImagem = selecionada ? (
    modo === 'gene' 
      ? `http://127.0.0.1:33857/plot_brain?gene=${selecionada.value}`
      : `http://127.0.0.1:33857/plot_ontology?geneset=${selecionada.value}`
  ) : '';

  return (
    <div className="container" style={{ fontFamily: 'sans-serif', maxWidth: '1000px', margin: '0 auto', padding: '20px' }}>
      <h1 style={{ textAlign: 'center', color: '#2c3e50' }}>BrainMap</h1>
      
      {/* Botões para alternar o modo */}
      <div style={{ display: 'flex', justifyContent: 'center', marginBottom: '20px', gap: '10px' }}>
        <button 
          onClick={() => setModo('ontologia')}
          style={{ 
            padding: '10px 20px', 
            backgroundColor: modo === 'ontologia' ? '#3498db' : '#ecf0f1',
            color: modo === 'ontologia' ? 'white' : '#2c3e50',
            border: 'none', borderRadius: '5px', cursor: 'pointer', fontWeight: 'bold'
          }}
        >
          Explorar Ontologias
        </button>
        <button 
          onClick={() => setModo('gene')}
          style={{ 
            padding: '10px 20px', 
            backgroundColor: modo === 'gene' ? '#3498db' : '#ecf0f1',
            color: modo === 'gene' ? 'white' : '#2c3e50',
            border: 'none', borderRadius: '5px', cursor: 'pointer', fontWeight: 'bold'
          }}
        >
          Explorar Genes
        </button>
      </div>
      
      {/* A nossa nova caixa de seleção turbinada */}
      <div style={{ marginBottom: '30px' }}>
        <label style={{ fontWeight: 'bold', display: 'block', marginBottom: '10px', color: '#34495e' }}>
          Pesquise ou escolha {modo === 'gene' ? 'um Gene' : 'uma Ontologia'}: 
        </label>
        
        <Select
          value={selecionada}
          onChange={(opcao) => {
            setSelecionada(opcao);
            setCarregandoImagem(true); // Quando troca a opção, liga o aviso de "Carregando"
          }}
          options={opcoes}
          isLoading={carregandoLista}
          isSearchable={true}
          placeholder="Digite para buscar..."
          noOptionsMessage={() => "Nenhum resultado encontrado"}
        />
      </div>

      {/* Área do Gráfico */}
      <div style={{ textAlign: 'center', minHeight: '500px' }}>
        
        {/* Mostra este texto enquanto o R está processando (pode trocar por um GIF ou Spinner depois) */}
        {carregandoImagem && (
          <h3 style={{ color: '#e67e22' }}>⏳ Desenhando o mapa cerebral... Por favor, aguarde.</h3>
        )}

        {selecionada && (
          <img 
            src={urlImagem} 
            alt={`Gráfico para ${selecionada.label}`} 
            // Oculta a imagem enquanto ela ainda está carregando para não ficar "piscando"
            style={{ 
              display: carregandoImagem ? 'none' : 'inline-block', 
              width: '100%', 
              maxWidth: '900px', 
              borderRadius: '8px',
              boxShadow: '0 4px 8px rgba(0,0,0,0.1)'
            }}
            // O SEGREDO: Quando a imagem termina de baixar, ele desliga o aviso de "Carregando"
            onLoad={() => setCarregandoImagem(false)} 
          />
        )}
      </div>
    </div>
  );
}

export default App;