import tkinter as tk
from tkinter import messagebox
from tkinter import scrolledtext
import random
import math
import json
import os
import numpy as np
from collections import Counter


def breaking_sequence(resultados, comprimento):
    resultados['Tamanhos'][comprimento] += 1
    resultados['Proteinas Interrompidas'] += 1

def simular_polimerizacao(resultados, num_proteinas, tamanho_maximo, pool, pesos_pool, tipos_aminos, pesos_tipos, aminoacidos_prot, aminoacidos_nao_prot, amines, probabilities, num=False, tam=False):
    if num and tam:
        count_completas = 0
        count_prot = 0
        resultados['Tamanhos'] = {i: 0 for i in range(0, tam + 1)}

    for _ in range(num_proteinas if not num else num):
        comprimento = 0
        possui_nao_proteinogenico = False
        todos_aminoacidos_L = True  # Verifica se todos os aminoácidos são L
        sequencia_hidropatica = True
        sequencia_boa = []
        sequencia_ruim = []

        while comprimento < (tamanho_maximo if not tam else tam):
            # Verificar sequencia hidropática
            hidrop = random.choices(tipos_aminos, weights=pesos_tipos, k=1)[0] if not hidropatica else None

            # Sortear um único aminoácido ou amina do pool
            sorteio = random.choices(pool, weights=pesos_pool, k=1)[0]

            if sorteio not in aminoacidos_prot:
                possui_nao_proteinogenico = True

            reac_lateral = False
            if sorteio in aminoacidos_prot or sorteio in aminoacidos_nao_prot:
                if (sorteio not in aminoacidos_prot or amino_prot_proportion[sorteio[-3:]][2] != hidrop) and comprimento > 1 and not hidropatica:
                    sequencia_hidropatica = False

                if todos_aminoacidos_L and racemico:
                    if sorteio.startswith('D') if sorteio in aminoacidos_prot else random.random() > 0.5:
                        todos_aminoacidos_L = False

                if (sorteio in amino_laterais['Padrão'] or sorteio in amino_laterais['Alternativos']) and laterais:
                    count_amidas = sum(aminoacido in grupos_reac['Amida'] for aminoacido in sequencia_ruim)
                    if random.randint(0, count_amidas + 2) > 1:
                        breaking_sequence(resultados, comprimento)
                        break
                    for grupo in amino_laterais['Padrão'][sorteio] if sorteio in aminoacidos_prot else \
                    amino_laterais['Alternativos'][sorteio]:
                        if grupo == 'Amida':  # Ligações Isopeptídica (com outro AA)
                            if random.randint(0, len(sequencia_ruim) + 2) > 1:
                                reac_lateral = True
                                break
                        count = sum(aminoacido in grupos_reac[grupo] for aminoacido in sequencia_ruim)
                        if random.randint(0, count + 2) > 1:
                            reac_lateral = True
                            break
                    if reac_lateral:
                        breaking_sequence(resultados, comprimento)
                        break

                sequencia_ruim.append(sorteio)
                sequencia_boa.append(sorteio) if not possui_nao_proteinogenico else None
                comprimento += 1

            elif sorteio in amines:
                if comprimento >= 1:  # Considera apenas se a proteína tem pelo menos 2 aminoácidos
                    breaking_sequence(resultados, comprimento)
                    break

        else:
            # Se o loop termina sem encontrar uma amina nem lig. lateral, considera que atingiu o tamanho máximo
            if num and tam:
                count_completas += 1
                if not possui_nao_proteinogenico:
                    count_prot += 1
            else:
                resultados['Proteinas Completas'] += 1
                resultados['Tamanhos'][comprimento] += 1
                if not possui_nao_proteinogenico:
                    resultados['Apenas Proteinogenicos'] += 1
                    resultados['sequencias'].append(sequencia_boa)
                else:
                    resultados['Com AA não Proteinogenicos'] += 1
                    resultados['sequencias_ruins'].append(sequencia_ruim)
                if sequencia_hidropatica:
                    resultados['Sequências Hidropáticas'] += 1
                if todos_aminoacidos_L:
                    resultados['Proteinas Homoquirais'] += 1  # Contando proteínas completas formadas por aminoácidos L
                if todos_aminoacidos_L and not possui_nao_proteinogenico:
                    resultados['Apenas Proteinogenicos Homoquirais'] += 1
                if todos_aminoacidos_L and not possui_nao_proteinogenico and sequencia_hidropatica:
                    resultados['Sequências Hidropáticas, Homoquirais e Proteinogênicas'] += 1

    if num and tam:
        return [count_completas / num, count_prot / num]
    return resultados

def simular_proteinas(percentuais, num_proteinas,
                      tamanho_maximo, aminoacidos_prot, pesos_prot, aminoacidos_nao_prot,
                      pesos_nao_prot, amines, pesos_amines, best_prot, total_homoquiral_chance):
    global prob_formar
    resultados = {
        'Tamanhos': {i: 0 for i in range(0, tamanho_maximo + 1)},
        'Proteinas Interrompidas': 0,
        'Proteinas Completas': 0,
        'Sequências Hidropáticas': 0,
        'Com AA não Proteinogenicos': 0,
        'Apenas Proteinogenicos': 0,
        'Proteinas Homoquirais': 0,
        'Apenas Proteinogenicos Homoquirais': 0,
        'Sequências Hidropáticas, Homoquirais e Proteinogênicas': 0,
        'sequencias': [],
        'sequencias_ruins': []
    }

    probabilities = {}
    # Verificando se todos os 20 aminoácidos estão presentes e têm a mesma abundância
    todos_presentes = len(base_amino_prot_proportion) == 20

    # Inicializando a probabilidade
    probabilities['Prob Sequência Útil'] = 0

    if todos_presentes:
        # Cálculo da probabilidade inicial
        probabilities['Prob Sequência Útil'] = 3.18 ** -tamanho_maximo * (
                    ((1 - (percentuais['percentual_aminas'])) ** tamanho_maximo) * (
                    total_homoquiral_chance ** tamanho_maximo) * (percentuais['percentual_proteinogenicos'] ** tamanho_maximo if percentuais['percentual_proteinogenicos'] > 0 else 1))

        # Encontre o máximo valor de abundância
        abundancias = list(base_amino_prot_proportion.values())
        max_abundancia = max(abundancias)

        # Verifica se algum aminoácido tem abundância igual a zero
        if any(abundancia == 0 for abundancia in abundancias):
            probabilities['Prob Sequência Útil'] = 0
        else:
            # Calcular fator de ajuste com base na proximidade a zero
            for abundancia in abundancias:
                if abundancia < max_abundancia:
                    proporcao = abundancia / max_abundancia
                    probabilities['Prob Sequência Útil'] *= proporcao

    # Criando o pool com todos os aminoácidos e aminas
    pool = list(aminoacidos_prot) + list(aminoacidos_nao_prot) + list(amines)
    pesos_pool = list(pesos_prot) + list(pesos_nao_prot) + list(pesos_amines)

    tipos_aminos = False
    pesos_tipos = False
    if not hidropatica:
        # Contar a frequência de cada item
        frequencia = Counter(values[2] for values in amino_prot_proportion.values())
        tipo_total = sum(frequencia.values())

        # Probabilidade de cada tipo de aminoácido
        tipo_probs = {tipo: freq / tipo_total for tipo, freq in frequencia.items()}

        # Separar os elementos e pesos para a seleção ponderada
        tipos_aminos = list(tipo_probs.keys())
        pesos_tipos = list(tipo_probs.values())

        pool_real = {
            'Apolar': 0,
            'Polar': 0,
            'Acid': 0,
            'Alcaline': 0
        }

        for i, aminoacido in enumerate(pool):
            if aminoacido[-3:] in amino_prot_proportion:
                tipo = amino_prot_proportion[aminoacido[-3:]][2]
                pool_real[tipo] += pesos_pool[i] / sum(pesos_pool)
        probabilities['Prob Sequências Hidropáticas'] = 0
        for key in pool_real:
            probabilities['Prob Sequências Hidropáticas'] += tipo_probs[key] * pool_real[key]
        probabilities['Prob Sequências Hidropáticas'] = probabilities['Prob Sequências Hidropáticas'] ** tamanho_maximo
        probabilities['Sequências Hidropáticas'] = probabilities['Prob Sequências Hidropáticas'] * num_proteinas
    else:
        probabilities['Prob Sequências Hidropáticas'] = 1

    resultados = simular_polimerizacao(resultados, num_proteinas, tamanho_maximo, pool, pesos_pool, tipos_aminos, pesos_tipos, aminoacidos_prot, aminoacidos_nao_prot, amines, probabilities)

    # Cálculo de porcentagens

    resultados['Tamanhos'] = {tamanho: round((valor / num_proteinas) * 100, 4)
    for tamanho, valor in resultados['Tamanhos'].items() if round((valor / num_proteinas) * 100, 4) != 0}
        #resultados['Tamanhos'][tamanho] = 0 if resultados['Tamanhos'][tamanho] == 0.0 else resultados['Tamanhos'][tamanho]

    resultados['% com AA não Proteinogenicos'] = (resultados['Com AA não Proteinogenicos'] / resultados[
        'Proteinas Completas']) * 100 if resultados['Proteinas Completas'] > 0 else 0
    resultados['% Apenas Proteinogenicos'] = 100 - resultados['% com AA não Proteinogenicos'] if \
        resultados['Proteinas Completas'] > 0 else 0

    # Cálculo do percentual de proteínas formadas apenas por aminoácidos proteinogênicos L
    resultados['% Apenas Proteinogenicos Homoquirais'] = (resultados['Apenas Proteinogenicos Homoquirais'] / resultados[
        'Proteinas Completas']) * 100 if resultados[ 'Proteinas Completas'] != 0 else 0

    # Inicializa variáveis para encontrar a melhor sequência
    melhor_numero = -float('inf')  # Valor inicial muito baixo
    melhor_sequencias = []

    if len(resultados['sequencias']) > 0:
        # Analisando todas as sequências
        for sequencia in resultados['sequencias'][:3]:
            total_aminoacidos = len(sequencia)
            total_comum = sum(sequencia.count(amino_acido) for amino_acido in best_prot)

            # Calculando o número restante de aminoácidos
            numero_restante = total_aminoacidos - total_comum

            # Verificando se o número restante é maior que o atual
            if numero_restante > melhor_numero:
                melhor_numero = numero_restante
                melhor_sequencias = [sequencia]  # Reinicia a lista com a nova melhor sequência
            elif numero_restante == melhor_numero:
                melhor_sequencias.append(sequencia)  # Adiciona à lista se é um empate
    else:
        # Analisando todas as sequências
        for sequencia in resultados['sequencias_ruins'][:3]:
            # total_aminoacidos = len(sequencia)

            # Contando a quantidade total de aminoácidos encontrados nas chaves
            total_amino_prot = sum(
                sequencia.count(amino_acido) for amino_acido in new_amino_prot_proportion.keys())

            # Verificando se o número restante é maior que o atual
            if total_amino_prot > melhor_numero:
                melhor_numero = total_amino_prot
                melhor_sequencias = [sequencia]  # Reinicia a lista com a nova melhor sequência
            elif total_amino_prot == melhor_numero:
                melhor_sequencias.append(sequencia)  # Adiciona à lista se é um empate

    # Calcular Projeções
    rota = "P("
    for amino_acido, entry in entries_prot.items():
        peso = "{:.1e}".format(float(entry.get().strip())).replace('.0e+00', '').replace('e+00', '').replace('.0e-0', '-').replace('e-0', '-')
        rota += peso + '/'
    rota = rota[:-1] + ")N("
    for amino_acido, entry in entries_nao_prot.items():
        peso = "{:.1e}".format(float(entry.get().strip())).replace('.0e+00', '').replace('e+00', '').replace('.0e-0', '-').replace('e-0', '-')
        rota += peso + '/'
    rota = rota[:-1] + ")A("
    for amine, entry in entries_amines.items():
        peso = "{:.1e}".format(float(entry.get().strip())).replace('.0e+00', '').replace('e+00', '').replace('.0e-0', '-').replace('e-0', '-')
        rota += peso + '/'
    rota = rota[:-1] + ")"

    caminho = f"{rota}{racemico}{laterais}"

    formou = resultados['Proteinas Completas'] / num_proteinas
    formou_prot = resultados['Apenas Proteinogenicos'] / num_proteinas

    if caminho not in prob_formar:
        prob_formar[caminho] = {
            tamanho_maximo: {
                'formacao': [],
                'apenas_prot': []
            }
        }

        if 10 not in prob_formar[caminho] or 11 not in prob_formar[caminho]:
            if 10 not in prob_formar[caminho]:
                prob_formar[caminho][10] = {
                    'formacao': [],
                    'apenas_prot': []
                }
            if 11 not in prob_formar[caminho]:
                prob_formar[caminho][11] = {
                    'formacao': [],
                    'apenas_prot': []
                }

            lengs, interrupt = resultados['Tamanhos'], resultados['Proteinas Interrompidas']
            for _ in range(100):
                sim_results = simular_polimerizacao(resultados, num_proteinas, tamanho_maximo, pool, pesos_pool, tipos_aminos, pesos_tipos, aminoacidos_prot, aminoacidos_nao_prot, amines, probabilities, 10000, 10)
                prob_formar[caminho][10]['formacao'].append(sim_results[0])
                prob_formar[caminho][10]['apenas_prot'].append(sim_results[1])
                sim_results = simular_polimerizacao(resultados, num_proteinas, tamanho_maximo, pool, pesos_pool, tipos_aminos, pesos_tipos, aminoacidos_prot, aminoacidos_nao_prot, amines, probabilities, 10000, 11)
                prob_formar[caminho][11]['formacao'].append(sim_results[0])
                prob_formar[caminho][11]['apenas_prot'].append(sim_results[1])
                print(f'calibração inicial: {_}% concluída')

            resultados['Tamanhos'] = lengs
            resultados['Proteinas Interrompidas'] = interrupt

    if tamanho_maximo not in prob_formar[caminho]:
        prob_formar[caminho][tamanho_maximo] = {
            'formacao': [],
            'apenas_prot': []
        }

    dif = (np.mean(prob_formar[caminho][11]['formacao']) / np.mean(prob_formar[caminho][10]['formacao']))
    dif_prot = (np.mean(prob_formar[caminho][11]['apenas_prot']) / np.mean(prob_formar[caminho][10]['apenas_prot']))

    if resultados['Proteinas Completas'] == 0:
        lista = []
        for key in prob_formar[caminho]:
            if key == tamanho_maximo:
                pass
            elif len(prob_formar[caminho][key]['formacao']) >= 10:
                lista.append(key)

        lista.sort(reverse=True)  # Ordena a lista em ordem crescente

        # Inicializa as variáveis para armazenar o menor intervalo e o par correspondente
        menor_intervalo = float('inf')
        par_mais_proximo = None

        # Itera sobre os pares consecutivos na lista
        for i in range(len(lista) - 1):
            intervalo = abs(lista[i + 1] - lista[i])
            if intervalo < menor_intervalo:
                menor_intervalo = intervalo
                par_mais_proximo = (lista[i + 1], lista[i])

        if tamanho_maximo > par_mais_proximo[1]:
            mean = np.mean(prob_formar[caminho][par_mais_proximo[1]]['formacao']) * (
                        dif ** (tamanho_maximo - par_mais_proximo[1]))
            mean_prot = np.mean(prob_formar[caminho][par_mais_proximo[1]]['apenas_prot']) * (
                        dif_prot ** (tamanho_maximo - par_mais_proximo[1]))
        else:
            mean = np.mean(prob_formar[caminho][par_mais_proximo[1]]['formacao']) / (
                        dif ** abs(tamanho_maximo - par_mais_proximo[1]))
            mean_prot = np.mean(prob_formar[caminho][par_mais_proximo[1]]['apenas_prot']) / (
                        dif_prot ** abs(tamanho_maximo - par_mais_proximo[1]))
        prob_formar[caminho][tamanho_maximo]['formacao'].append(mean)
        prob_formar[caminho][tamanho_maximo]['apenas_prot'].append(mean_prot)
    else:
        prob_formar[caminho][tamanho_maximo]['formacao'].append(formou)
        prob_formar[caminho][tamanho_maximo]['apenas_prot'].append(formou_prot)
        if len(prob_formar[caminho][tamanho_maximo]['formacao']) > 9:
            mean = np.mean(prob_formar[caminho][tamanho_maximo]['formacao'])
            mean_prot = np.mean(prob_formar[caminho][tamanho_maximo]['apenas_prot'])
        else:
            mean11 = (np.mean(prob_formar[caminho][11]['formacao']) * (dif ** abs(tamanho_maximo - 11))) * num_proteinas
            meanprot11 = (np.mean(prob_formar[caminho][11]['apenas_prot']) * (dif ** abs(tamanho_maximo - 11))) * num_proteinas
            if mean11/2 > resultados['Proteinas Completas'] > mean11*2:
                mean = formou
            else:
                mean = np.mean(prob_formar[caminho][11]['formacao']) * (dif ** abs(tamanho_maximo - 11))
            if meanprot11/2 > resultados['Apenas Proteinogenicos'] > meanprot11*2:
                mean_prot = formou
            else:
                mean_prot = np.mean(prob_formar[caminho][11]['apenas_prot']) * (dif_prot ** abs(tamanho_maximo - 11))

    # Salva o resultado da melhor sequência em 'resultados'
    resultados['melhor_sequencia'] = melhor_sequencias

    with open('prob_formar.json', 'w') as file:
        json.dump(prob_formar, file, indent=4)

    probabilities['Prob de Formar'] = mean if mean else formou
    probabilities['Proteinas Completas'] = probabilities['Prob de Formar'] * num_proteinas
    probabilities['Sequências Hidropáticas'] = probabilities['Proteinas Completas'] if hidropatica else probabilities['Sequências Hidropáticas']
    probabilities['Prob Apenas Proteinogenicos'] = mean_prot if mean_prot else ((1 - (percentuais['percentual_aminas'])) ** tamanho_maximo) * (
            percentuais['percentual_proteinogenicos'] ** tamanho_maximo)
    probabilities['Apenas Proteinogenicos'] = probabilities['Prob Apenas Proteinogenicos'] * num_proteinas
    probabilities['% Apenas Proteinogenicos'] = (probabilities['Apenas Proteinogenicos'] / probabilities['Proteinas Completas']) * 100
    print(percentuais['glycine'])
    probabilities['Prob Proteinas Homoquirais'] = probabilities['Prob de Formar'] * (1 + percentuais['glycine']) * (total_homoquiral_chance ** tamanho_maximo) if racemico else 1
    probabilities['Proteinas Homoquirais'] = probabilities['Prob Proteinas Homoquirais'] * num_proteinas * (1 + percentuais['glycine'])
    probabilities['Prob Apenas Proteinogenicos Homoquirais'] = ((total_homoquiral_chance ** tamanho_maximo) * (1 + percentuais['glycine']) if racemico else 1) * probabilities['Prob Apenas Proteinogenicos']
    probabilities['Apenas Proteinogenicos Homoquirais'] = probabilities['Prob Apenas Proteinogenicos Homoquirais'] * num_proteinas
    probabilities['% Apenas Proteinogenicos Homoquirais'] = probabilities['% Apenas Proteinogenicos'] *  ((total_homoquiral_chance ** tamanho_maximo) if racemico else 1)
    probabilities['Prob Sequências Hidropáticas, Homoquirais e Proteinogênicas'] = (mean_prot if mean_prot else ((1 - (percentuais['percentual_aminas'])) ** tamanho_maximo) * (
            percentuais['percentual_proteinogenicos'] ** tamanho_maximo)) * probabilities['Prob Sequências Hidropáticas'] * (((total_homoquiral_chance ** tamanho_maximo) * (1 + percentuais['glycine'])) if racemico else 1)
    probabilities['Sequências Hidropáticas, Homoquirais e Proteinogênicas'] = probabilities['Prob Sequências Hidropáticas, Homoquirais e Proteinogênicas'] * num_proteinas
    probabilities['Prob Sequência Útil'] = (3.18 ** -tamanho_maximo) * probabilities['Prob Sequências Hidropáticas, Homoquirais e Proteinogênicas'] if todos_presentes else 0
    probabilities['Sequência Útil'] = probabilities['Prob Sequência Útil'] * num_proteinas

    for key in probabilities.keys():
        if key.startswith('Prob'):
            print(f"{key}: {probabilities[key]}")
    print('--------')
    return resultados, melhor_numero, probabilities


def calcular_proporcoes():
    global percentuais, base_amino_prot_proportion, new_amino_prot_proportion
    # Obtem valores dos campos de entrada e ajusta os dicionários

    #try:
    new_amino_prot_proportion = {}
    new_amino_nao_prot_proportion = {}
    new_amines_proportion = {}
    base_amino_prot_proportion = {}

    # Populando new_amino_prot_proportion com as quantidades de isômeros
    if racemico:
        for amino_acido, entry in entries_prot.items():
            value = float(entry.get().strip())  # Assume que o segundo elemento é o valor
            if value and float(value) > 0:
                if amino_acido == 'GLY':
                    new_amino_prot_proportion[amino_acido] = value  # Deixa como está
                    base_amino_prot_proportion[amino_acido] = value
                else:
                    # Calcula as proporções para L e D
                    if amino_acido == 'ILE' or amino_acido == 'THR':  # Para ILE ou THR
                        new_amino_prot_proportion[f'L-{amino_acido}'] = value / 4
                        new_amino_prot_proportion[f'D-{amino_acido}'] = (value * 3) / 4
                        base_amino_prot_proportion[amino_acido] = value
                    else:  # Para outros aminoácidos com 2 isômeros
                        new_amino_prot_proportion[f'L-{amino_acido}'] = value / 2
                        new_amino_prot_proportion[f'D-{amino_acido}'] = value / 2
                        base_amino_prot_proportion[amino_acido] = value
    else:
        for amino_acido, entry in entries_prot.items():
            value = entry.get().strip()  # Remove espaços em branco
            if value and float(value) > 0:
                new_amino_prot_proportion[amino_acido] = float(value)  # Se estiver vazio, considera 0
        base_amino_prot_proportion = new_amino_prot_proportion
    # Coleta os valores dos aminoácidos não proteogênicos
    for amino_acido, entry in entries_nao_prot.items():
        value = entry.get().strip()  # Remove espaços em branco
        new_amino_nao_prot_proportion[amino_acido] = float(value) if value else 0.0  # Se estiver vazio, considera 0

    # Coleta os valores das aminas
    for amino_acido, entry in entries_amines.items():
        value = entry.get().strip()  # Remove espaços em branco
        new_amines_proportion[amino_acido] = float(value) if value else 0.0  # Se estiver vazio, considera 0

    num_proteinas = int(entry_num_proteinas.get())
    tamanho_maximo = int(entry_tamanho_maximo.get())

    # Calculando totais
    total_prot = sum(new_amino_prot_proportion.values())
    total_nao_prot = sum(new_amino_nao_prot_proportion.values())
    total_amines = sum(new_amines_proportion.values())

    # Total combinado
    total_geral = total_prot + total_nao_prot + total_amines
    total_aminos = total_prot + total_nao_prot

    # Criando um dicionário para os percentuais de cada tipo de ligação
    percentuais = {}

    percentuais['percentual_proteinogenicos'] = total_prot / total_geral
    percentuais['percentual_nao_proteinogenicos'] = total_nao_prot / total_geral
    percentuais['percentual_aminas'] = total_amines / total_geral
    percentuais['glycine'] = new_amino_prot_proportion['GLY'] / total_aminos
    percentuais['amida_prot'] = sum(new_amino_prot_proportion[amino] for amino in lig_laterais['Amida'][:9] if amino in new_amino_prot_proportion) / total_aminos
    percentuais['amida_nao_prot'] = sum(new_amino_nao_prot_proportion[amino] for amino in lig_laterais['Amida'][9:] if amino in new_amino_nao_prot_proportion) / total_aminos
    percentuais['react_prot'] = 0
    percentuais['react_nao_prot'] = 0
    for tipo in lig_laterais:
        if tipo == 'Hidroxila':
            percentuais['react_prot'] += sum(new_amino_prot_proportion[amino] for amino in lig_laterais[tipo][:9] if amino in new_amino_prot_proportion) / total_aminos
            percentuais['react_nao_prot'] += sum(new_amino_nao_prot_proportion[amino] for amino in lig_laterais[tipo][9:] if amino in new_amino_nao_prot_proportion) / total_aminos
        elif tipo == 'Fenol':
            percentuais['react_prot'] += sum(new_amino_prot_proportion[amino] for amino in lig_laterais[tipo] if amino in new_amino_prot_proportion) / total_aminos
        elif tipo == 'Tiol':
            percentuais['react_prot'] += sum(new_amino_prot_proportion[amino] for amino in lig_laterais[tipo][:6] if amino in new_amino_prot_proportion) / total_aminos
            percentuais['react_nao_prot'] += sum(new_amino_nao_prot_proportion[amino] for amino in lig_laterais[tipo][6:] if amino in new_amino_nao_prot_proportion) / total_aminos
    for prob in percentuais:
        print(f"{prob}: {percentuais[prob]}")

    for tipo_ligacao, amino_acidos in lig_laterais.items():
        percentual_total = 0

        for amino_acido in amino_acidos:
            # Verifica se o aminoácido está nos dicionários e calcula a proporção
            if amino_acido in new_amino_prot_proportion:
                percentual_total += new_amino_prot_proportion[amino_acido] / total_geral
            elif amino_acido in new_amino_nao_prot_proportion:
                percentual_total += new_amino_nao_prot_proportion[amino_acido] / total_geral
            elif amino_acido in new_amines_proportion:
                percentual_total += new_amines_proportion[amino_acido] / total_geral

                # Armazenando o percentual total para cada tipo de ligação
        percentuais[f'percentual_{tipo_ligacao}'] = percentual_total
    # Percentuais de aminoácidos que têm ligações laterais
    percentuais_com_ligacao = sum([percentuais[f'percentual_{tipo}'] for tipo in lig_laterais])

    # Calculando percentual sem ligações
    percentuais[f'percentual_sem_lig'] = (percentuais['percentual_proteinogenicos'] + percentuais['percentual_nao_proteinogenicos'] + percentuais['percentual_aminas'] - percentuais_com_ligacao)

    aminoacidos_prot, pesos_prot = zip(
        *[(amino_acido, peso) for amino_acido, peso in new_amino_prot_proportion.items() if peso != 0])

    if total_nao_prot != 0:
        aminoacidos_nao_prot, pesos_nao_prot = zip(
            *[(amino_acido, peso) for amino_acido, peso in new_amino_nao_prot_proportion.items() if peso != 0])
    else:
        aminoacidos_nao_prot, pesos_nao_prot = [], []

    if total_amines != 0:
        amines, pesos_amines = zip(
            *[(amine, peso) for amine, peso in new_amines_proportion.items() if peso != 0])
    else:
        amines, pesos_amines = [], []

    # Inicializando variáveis
    total_weight = 0
    total_homoquiral_weight = 0

    # Cálculo da probabilidade de ser homoquiral
    for isomer, abundance in new_amino_prot_proportion.items():
        total_weight += abundance  # Soma total de abundâncias
        if isomer.startswith('L-') or isomer == 'GLY':
            total_homoquiral_weight += abundance  # Soma apenas as abundâncias L

    # A porcentagem de homoquiralidade é a proporção de isômeros L no total
    if total_weight > 0:
        percentual_prot_homoquiral = (total_homoquiral_weight / total_weight) * 100
    else:
        percentual_prot_homoquiral = 0

    total_homoquiral_chance = ((percentual_prot_homoquiral * percentuais['percentual_proteinogenicos']) + (
                50 * percentuais['percentual_nao_proteinogenicos'])) / (
                                          percentuais['percentual_proteinogenicos'] + percentuais['percentual_nao_proteinogenicos']) / 100 if racemico else 1

    # Encontrar o maior valor em ambos os dicionários
    max_value_prot = max(new_amino_prot_proportion.values())
    max_value_nao_prot = max(new_amino_nao_prot_proportion.values())

    # Determinar o maior valor entre os dois
    max_value = max(max_value_prot, max_value_nao_prot)

    # Calcular 10% do maior valor
    threshold = 0.2 * max_value

    # Criar uma lista com as chaves cujo valores são iguais ou maiores que o limite
    best_prot = [
        key for key, value in new_amino_prot_proportion.items() if value >= threshold
    ]
    best_nao_prot = [
        key for key, value in new_amino_nao_prot_proportion.items() if value >= threshold
    ]

    # Executa a simulação
    resultados, numero_restante, probabilities = simular_proteinas(percentuais, num_proteinas, tamanho_maximo,
                                                                   aminoacidos_prot, pesos_prot, aminoacidos_nao_prot,
                                                                   pesos_nao_prot, amines, pesos_amines, best_prot,
                                                                   total_homoquiral_chance)

    # Exibindo resultados na caixa de texto
    result_box.delete(1.0, tk.END)
    result_box.insert(tk.END, "RESULTADOS DA SIMULAÇÃO:\n\n", 'title')

    # Mapeamento de formatação
    formato_aminoacidos = {}

    # Adicionar formatação para aminoácidos proteinogênicos
    for aminoacido in new_amino_prot_proportion:
        if aminoacido not in best_prot:
            if aminoacido.startswith('D'):
                formato_aminoacidos[aminoacido] = "blue under"
            else:
                formato_aminoacidos[aminoacido] = "blue"
        else:
            if aminoacido.startswith('D'):
                formato_aminoacidos[aminoacido] = "under"
            else:
                formato_aminoacidos[aminoacido] = "normal"  # ou None

    # Adicionar formatação para aminoácidos não proteinogênicos
    for aminoacido in new_amino_nao_prot_proportion:
        formato_aminoacidos[aminoacido] = "red"

    for i, key in enumerate(resultados.keys()):
        if i == 11:
            result_box.insert(tk.END, '\n')
            result_box.insert(tk.END, f"{key} ", "bold")  # Inserir a key em negrito
            result_box.insert(tk.END, f"{round(resultados[key], 2)}%\n", 'normal')  # Inserir o resultado normal
        elif key == 'sequencias':
            pass
        elif key == 'sequencias_ruins':
            pass
        elif key == 'melhor_sequencia':
            if len(resultados[key]) > 1:
                result_box.insert(tk.END, '\nMelhores Sequências: ', 'bold')
                result_box.insert(tk.END, str(numero_restante) + '\n', 'normal')
            elif len(resultados[key]) == 1:
                result_box.insert(tk.END, '\nMelhor Sequência: ', 'bold')
                result_box.insert(tk.END, str(numero_restante) + '\n', 'normal')
            else:
                pass
            for j, sublist in enumerate(resultados[key]):
                result_box.insert(tk.END, str(j + 1) + '. {', 'normal')
                for i, amino_acido in enumerate(sublist):
                    if i > 0:
                        result_box.insert(tk.END, ", ", 'normal')
                    tag = formato_aminoacidos.get(amino_acido, None)  # Pega a formatação, se houver
                    if tag == "blue":
                        result_box.insert(tk.END, f"{amino_acido}", "blue")
                    elif tag == 'blue under':
                        result_box.insert(tk.END, f"{amino_acido}", "blue under")
                    elif tag == 'under':
                        result_box.insert(tk.END, f"{amino_acido}", "under")
                    elif tag == "red":
                        result_box.insert(tk.END, f"{amino_acido}", "red")
                    else:
                        result_box.insert(tk.END, f"{amino_acido}", 'normal')  # Formato normal

                result_box.insert(tk.END, "}\n\n")  # Nova linha após cada sublista
        else:
            result_box.insert(tk.END, f"{key} ", "bold")  # Inserir a key em negrito
            if key in list(probabilities.keys()):
                result_box.insert(tk.END, f"({round(probabilities[key], 4)})", 'red')
            if isinstance(resultados[key], float):
                result_box.insert(tk.END, f": {round(resultados[key], 2)}%\n", 'normal')  # Inserir o resultado normal
            else:
                result_box.insert(tk.END, f": {resultados[key]}\n", 'normal')  # Inserir o resultado normal

    # Adiciona Estatísticas comparativas:
    probabilidade = probabilities['Prob Sequência Útil']

    if probabilidade == 0:
        result_box.insert(tk.END, "-------------------------\nÉ impossível formar qualquer proteína útil.\n", 'normal')
        return

        # Converte a probabilidade para notação científica
    ordem_magnitude = math.floor(math.log10(probabilidade))
    valor_base = probabilidade / (10 ** ordem_magnitude)


    # Exibe os resultados em notação científica
    result_box.insert(tk.END, "-------------------------\nProbabilidade de uma Sequência Útil: ", 'bold')
    result_box.insert(tk.END,f"{valor_base:.2f} x 10", 'normal')
    result_box.insert(tk.END, f"{int(ordem_magnitude)}\n", 'superscript')

    #
    nome_prot = 'proteína' if tamanho_maximo > 49 else ('oligopeptídeo' if tamanho_maximo > 29 else 'oligopeptídeo curto')

    insert_weight(result_box, ordem_magnitude, valor_base, nome_prot)

    #except ValueError:
    #    messagebox.showerror("Erro", "Por favor, insira valores válidos.")

def insert_weight(result_box, ordem_magnitude, valor_base, nome_prot):
    if ordem_magnitude <= -80:
        if ordem_magnitude <= -320:
            grau = (-ordem_magnitude - 319) // 3
            expo = 319 + (grau * 3)
        elif ordem_magnitude <= -240:
            grau = (-ordem_magnitude - 239) // 3
            expo = 239 + (grau * 3)
        elif ordem_magnitude <= -160:
            grau = (-ordem_magnitude - 159) // 3
            expo = 159 + (grau * 3)
        else:
            grau = (-ordem_magnitude - 79) // 3
            expo = 79 + (grau * 3)
        peso_inutil = (1 / valor_base) * (10 ** (abs(ordem_magnitude) - expo))
        if peso_inutil < 1:
            peso_inutil *= 1000
            grau -= 1
        peso_str = f"{peso_inutil:.2f} {nomemclaturas[grau]}"
        comparacao = 'o número de átomos do Universo'
    elif ordem_magnitude <= -73:
        grau = (-ordem_magnitude - 72) // 3
        expo = 72 + (grau * 3)
        peso_inutil = (1 / valor_base) * (10 ** (abs(ordem_magnitude) - expo))
        if peso_inutil < 1:
            peso_inutil *= 1000
            grau -= 1
        peso_str = f"{peso_inutil/1.673864283083928:.2f} {nomemclaturas[grau]}"
        comparacao = 'a massa do Universo'
    elif ordem_magnitude <= -59:
        grau = (-ordem_magnitude - 58) // 3
        expo = 58 + (grau * 3)
        peso_inutil = (1 / valor_base) * (10 ** (abs(ordem_magnitude) - expo))
        if peso_inutil < 1:
            peso_inutil *= 1000
            grau -= 1
        peso_str = f"{peso_inutil/1.005570862578686:.2f} {nomemclaturas[grau]}"
        comparacao = 'a massa da Via Láctea'
    elif ordem_magnitude <= -50:
        grau = (-ordem_magnitude - 49) // 3
        expo = 49 + (grau * 3)
        peso_inutil = (1 / valor_base) * (10 ** (abs(ordem_magnitude) - expo))
        if peso_inutil < 1:
            peso_inutil *= 1000
            grau -= 1
        peso_str = f"{peso_inutil/5.02785431289343:.2f} {nomemclaturas[grau]}"
        comparacao = 'a massa do Sol'
    elif ordem_magnitude <= -47:
        peso_inutil = (1 / valor_base) * (10 ** (abs(ordem_magnitude) - 46))
        if peso_inutil < 1:
            peso_inutil *= 1000
        peso_str = f"{peso_inutil/5.268703898840885:.2f}"
        comparacao = 'a massa de Jupter'
    elif ordem_magnitude <= -44:
        peso_inutil = (1 / valor_base) * (10 ** (abs(ordem_magnitude) - 43))
        if peso_inutil < 1:
            peso_inutil *= 1000
        peso_str = f"{peso_inutil / 1.675041876046901:.2f}"
        comparacao = 'a massa da Terra'
    elif ordem_magnitude <= -42:
        peso_inutil = (1 / valor_base) * (10 ** (abs(ordem_magnitude) - 41))
        if peso_inutil < 1:
            peso_inutil *= 1000
        peso_str = f"{peso_inutil / 1.37438953472:.2f}"
        comparacao = 'a massa da Lua'
    elif ordem_magnitude <= -40:
        peso_inutil = (1 / valor_base) * (10 ** (abs(ordem_magnitude) - 39))
        if peso_inutil < 1:
            peso_inutil *= 1000
        peso_str = f"{peso_inutil / 1.073741824:.2f}"
        comparacao = 'a massa do Everest'
    elif ordem_magnitude <= -34:
        grau = (-ordem_magnitude - 33) // 3
        expo = 33 + (grau * 3)
        peso_inutil = (1 / valor_base) * (10 ** (abs(ordem_magnitude) - expo))
        if peso_inutil < 1:
            peso_inutil *= 1000
            grau -= 1
        peso_str = f"{peso_inutil / 4.19:.2f} {nomemclaturas[grau]}"
        comparacao = 'a massa de todos os seres humanos da Terra'
    elif ordem_magnitude <= -26:
        grau = (-ordem_magnitude - 25) // 3
        expo = 25 + (grau * 3)
        peso_inutil = (1 / valor_base) * (10 ** (abs(ordem_magnitude) - expo))
        if peso_inutil < 1:
            peso_inutil *= 1000
            grau -= 1
        peso_str = f"{peso_inutil:.2f} {nomemclaturas[grau]}"
        comparacao = 'Toneladas'
    elif ordem_magnitude <= -23:
        peso_inutil = (1 / valor_base) * (10 ** (abs(ordem_magnitude) - 5))
        if peso_inutil < 1:
            peso_inutil *= 1000
        peso_str = f"{(peso_inutil/6.02)*5.5:.2f}"
        comparacao = 'Kg'
    elif ordem_magnitude <= -22:
        peso_inutil = (1 / valor_base) * (10 ** (abs(ordem_magnitude) - 21))
        if peso_inutil < 1:
            peso_inutil *= 1000
        peso_str = f"{peso_inutil / 1.37438953472:.2f}"
        comparacao = 'a massa da Lua'
    elif ordem_magnitude <= -20:
        peso_inutil = (1 / valor_base) * (10 ** (abs(ordem_magnitude) - 19))
        if peso_inutil < 1:
            peso_inutil *= 1000
        peso_str = f"{peso_inutil / 1.073741824:.2f}"
        comparacao = 'a massa do Everest'
    elif ordem_magnitude <= -14:
        grau = (-ordem_magnitude - 13) // 3
        expo = 13 + (grau * 3)
        peso_inutil = (1 / valor_base) * (10 ** (abs(ordem_magnitude) - expo))
        if peso_inutil < 1:
            peso_inutil *= 1000
            grau -= 1
        peso_str = f"{peso_inutil / 4.19:.2f} {nomemclaturas[grau]}"
        comparacao = 'a massa de todos os seres humanos da Terra'
    elif ordem_magnitude <= -6:
        grau = (-ordem_magnitude - 5) // 3
        expo = 5 + (grau * 3)
        peso_inutil = (1 / valor_base) * (10 ** (abs(ordem_magnitude) - expo))
        if peso_inutil < 1:
            peso_inutil *= 1000
            grau -= 1
        peso_str = f"{peso_inutil:.2f} {nomemclaturas[grau]}"
        comparacao = 'Toneladas'
    elif ordem_magnitude <= -3:
        peso_inutil = (1 / valor_base) * (10 ** (abs(ordem_magnitude) - 5))
        if peso_inutil < 1:
            peso_inutil *= 1000
        peso_str = f"{peso_inutil:.2f}"
        comparacao = 'Kg'

    if ordem_magnitude <= -23:
        comp_prot = 'molécula'
    else:
        comp_prot = 'grama'

    if ordem_magnitude > -160:# Formatando partes da frase
        result_box.insert(tk.END, "Há aproximadamente ", 'normal')
        result_box.insert(tk.END, f"{peso_str} ", 'bold')
        result_box.insert(tk.END, "vezes ", 'normal') if comparacao != "Toneladas" and comparacao != 'Kg' else None
        result_box.insert(tk.END, f"{comparacao} ", 'bold')
        result_box.insert(tk.END, "de proteínas deletérias para cada ", 'normal')
        result_box.insert(tk.END, f"{comp_prot}", 'bold')
        result_box.insert(tk.END, f" de {nome_prot} útil.\n", 'normal')
    elif ordem_magnitude > -240:# Formatando partes da frase
        result_box.insert(tk.END, "Teria que transformar ", 'normal')
        result_box.insert(tk.END, "cada átomo do universo ", 'bold')
        result_box.insert(tk.END, "em proteínas deletérias ", 'normal')
        result_box.insert(tk.END, f"{peso_str} ", 'bold')
        result_box.insert(tk.END, "vezes, ", 'normal')
        result_box.insert(tk.END, "para cada ", 'normal')
        result_box.insert(tk.END, f"{comp_prot}", 'bold')
        result_box.insert(tk.END, f" de {nome_prot} útil.\n", 'normal')
    elif ordem_magnitude > -320:
        result_box.insert(tk.END, "Teria que transformar ", 'normal')
        result_box.insert(tk.END, "cada átomo do universo ", 'bold')
        result_box.insert(tk.END, "em outro universo e então transformar ", "normal")
        result_box.insert(tk.END, "cada átomo de cada um desses universos ", "bold")
        result_box.insert(tk.END, "em proteínas deletérias ", 'normal')
        result_box.insert(tk.END, f"{peso_str} ", 'bold')
        result_box.insert(tk.END, "vezes, ", 'normal')
        result_box.insert(tk.END, "para cada ", 'normal')
        result_box.insert(tk.END, f"{comp_prot}", 'bold')
        result_box.insert(tk.END, f" de {nome_prot} útil.\n", 'normal')
    else:
        result_box.insert(tk.END, "Teria que transformar ", 'normal')
        result_box.insert(tk.END, "cada átomo do universo ", 'bold')
        result_box.insert(tk.END, "em outro universo e então transformar ", "normal")
        result_box.insert(tk.END, "cada átomo de cada um desses outros universos ", 'bold')
        result_box.insert(tk.END, "em ainda outro universo, ", 'normal')
        result_box.insert(tk.END, "e repetir esse processo mais ", 'normal') if (ordem_magnitude // 80) - 3 > 1 else None
        result_box.insert(tk.END, f"{(ordem_magnitude // 80) - 3} vez", 'bold') if (ordem_magnitude // 80) - 3 > 1 else None
        result_box.insert(tk.END, "para só então transformar  ", 'normal')
        result_box.insert(tk.END, "todos esse átomos de todos esses universos", 'bold')
        result_box.insert(tk.END, "proteínas deletérias", 'normal')
        result_box.insert(tk.END, f"{peso_str} ", 'bold')
        result_box.insert(tk.END, "vezes, ", 'normal')
        result_box.insert(tk.END,
                          "e então transformar cada átomo de todos esses universos em proteínas deletérias para cada ",
                          'normal')
        result_box.insert(tk.END, f"{comp_prot}", 'bold')
        result_box.insert(tk.END, f" de {nome_prot} útil.\n", 'normal')

# Lista de nomenclaturas de grandezas
nomemclaturas = ['', 'mil ', 'milhões de ', 'bilhões de', 'trilhões de', 'quadrilhões de', 'quintilhões de',
                     'sextilhões de', 'septilhões de', 'octilhões de', 'nonilhões de', 'decilhõões de',
                     'undecilhões de', 'dodecilhões de', 'tridecilhões de', 'quadricilhões de', 'quinticilhões de',
                     'sexicilhões de ', 'septicilhões de ', 'octicilhões de ', 'nonicilhões de ', 'eicocilliões de ',
                     'uneicocilliões de ', 'dueicocilliões de ', 'trieicocilliões de ', 'quadricocilliões de ', 'quinticocilliões de ']

# Dicionário dos aminoácidos proteogênicos com nome e proporção
amino_prot_proportion = {
    'GLY': ['Glycine', 1, 'Apolar'],
    'ALA': ['Alanine', 1.1, 'Apolar'],
    'SER': ['Serine', 0.05, 'Polar'],
    'THR': ['Threonine', 0.02, 'Polar'],
    'VAL': ['Valine', 0.6, 'Apolar'],
    'ILE': ['Isoleucine', 0.03, 'Apolar'],
    'LEU': ['Leucine', 0.009, 'Apolar'],
    'MET': ['Methionine', 0.002, 'Apolar'],
    'PHE': ['Phenylalanine', 0, 'Apolar'],
    'TYR': ['Tyrosine', 0, 'Polar'],
    'TRP': ['Tryptophan', 0, 'Apolar'],
    'ASP': ['Aspartic Acid', 0.008, 'Acid'],
    'GLU': ['Glutamic Acid', 0.02, 'Acid'],
    'ASN': ['Asparagine', 0, 'Polar'],
    'GLN': ['Glutamine', 0, 'Polar'],
    'CYS': ['Cysteine', 0, 'Polar'],
    'LYS': ['Lysine', 0, 'Alcaline'],
    'ARG': ['Arginine', 0, 'Alcaline'],
    'HIS': ['Histidine', 0, 'Alcaline'],
    'PRO': ['Proline', 0, 'Apolar']
}

# Dicionário dos aminoácidos não proteinogênicos com nome e proporção
amino_nao_prot_proportion = {
    'B-ALA': ('β-Alanine', 0.9),
    'ISS': ('Isoserine', 0.003),
    'HCA': ('Homocysteic Acid', 0.0009),
    'B-ABA': ('β-Aminobutyric Acid', 0.7),
    'A-ABA': ('α-Aminobutyric Acid', 0.4),
    'A-AIB': ('α-Aminoisobutyric Acid', 0.1),
    'G-ABA': ('γ-Aminobutyric Acid', 0.1),
    'B-AIB': ('β-Aminoisobutyric Acid', 0.04),
    'SMC': ('S-Methylcysteine', 0.0002),
    'ISV': ('Isovaline', 0.095),
    'MSO': ('Methionine Sulfoxide', 0.00009),
    'MSF': ('Methionine Sulfone', 0.00002),
    'ETH': ('Ethionine', 0.003)
}


# Dicionário das aminas com nome e proporção
amines_proportion = {
    'MAM': ('Methylamine', 0.02),
    'EAM': ('Ethylamine', 0.05),
    'CYM': ('Cysteamine', 0.05),
    'ETHA': ('Ethanolamine', 0.04)
}

amino_laterais = {
    'Padrão': {
        'SER': ['Hidroxila'],
        'THR': ['Hidroxila'],
        'TYR': ['Hidroxila', 'Fenol'],
        'CYS': ['Tiol'],
        'MET': ['Tiol'],
        'LYS': ['Amida'],
        'ASP': ['Amida'],
        'GLU': ['Amida'],
        #'ARG': ['Guanidino'],
        #'HIS': ['Imidazol'],
        #'TRP': ['Indol'],
        'L-SER': ['Hidroxila'],
        'L-THR': ['Hidroxila'],
        'L-TYR': ['Hidroxila', 'Fenol'],
        'L-CYS': ['Tiol'],
        'L-MET': ['Tiol'],
        'L-LYS': ['Amida'],
        'L-ASP': ['Amida'],
        'L-GLU': ['Amida'],
        #'L-ARG': ['Guanidino'],
        #'L-HIS': ['Imidazol'],
        #'L-TRP': ['Indol'],
        'D-SER': ['Hidroxila'],
        'D-THR': ['Hidroxila'],
        'D-TYR': ['Hidroxila', 'Fenol'],
        'D-CYS': ['Tiol'],
        'D-MET': ['Tiol'],
        'D-LYS': ['Amida'],
        'D-ASP': ['Amida'],
        'D-GLU': ['Amida'],
        #'D-ARG': ['Guanidino'],
        #'D-HIS': ['Imidazol'],
        #'D-TRP': ['Indol'],
    },
    'Alternativos': {
        'ISO': ['Hidroxila', 'Amida'],
        'B-ABA': ['Amida'],
        'A-ABA': ['Amida'],
        'A-AIB': ['Amida'],
        'G-ABA': ['Amida'],
        'B-AIB': ['Amida'],
        'SMC': ['Tiol'],
        'HCA': ['Tiol'],
    }
}

# Criando o dicionário oposto
lig_laterais = {
    'Hidroxila': ['SER', 'THR', 'TYR', 'L-SER', 'L-THR', 'L-TYR', 'D-SER', 'D-THR', 'D-TYR', 'ISO'],
    'Fenol': ['TYR', 'L-TYR', 'D-TYR'],
    'Tiol': ['CYS', 'MET', 'L-CYS', 'L-MET', 'D-CYS', 'D-MET', 'SMC', 'HCA'],
    'Amida': ['LYS', 'ASP', 'GLU', 'L-LYS', 'L-ASP', 'L-GLU', 'D-LYS', 'D-ASP', 'D-GLU', 'ISO', 'B-ABA', 'A-ABA', 'A-AIB', 'G-ABA', 'B-AIB']
}

grupos_reac = {
        'Hidroxila': ['SER', 'THR', 'TYR', 'ASP', 'GLU', 'L-SER', 'L-THR', 'L-TYR', 'L-ASP', 'L-GLU', 'D-SER', 'D-THR', 'D-TYR', 'D-ASP', 'D-GLU', 'ISV'],  # Esterificação (com outro Hidroxila ou carboxila), Fosforilação (com -HPO4) e Glicosilação (com açúcares)
        'Tiol': ['CYS', 'MET', 'ASP', 'GLU', 'L-CYS', 'L-MET', 'L-ASP', 'L-GLU', 'D-CYS', 'D-MET', 'D-ASP', 'D-GLU', 'ISV'],  # Ligações Dissulfeto -S-S-, Tioésteres ou Tioéteres
        'Amida': ['SER', 'THR', 'TYR', 'L-SER', 'L-THR', 'L-TYR', 'D-SER', 'D-THR', 'D-TYR', 'ISV'],  # Esterificação (com outro Hidroxila)
        'Fenol': ['SER', 'THR', 'TYR', 'L-SER', 'L-THR', 'L-TYR', 'D-SER', 'D-THR', 'D-TYR']  # Esterificação (com outro Hidroxila)]
    }

# Carregar prob_formar do arquivo no início
if os.path.exists('prob_formar.json'):
    with open('prob_formar.json', 'r') as file:
        prob_formar = json.load(file)
    # Converter chaves internas para int
    for outer_key, inner_dict in prob_formar.items():
        if isinstance(inner_dict, dict):  # Verifica se o valor é um dicionário
            prob_formar[outer_key] = {int(inner_key): value for inner_key, value in inner_dict.items()}

else:
    prob_formar = {}  # Cria um dicionário vazio se o arquivo não existir

def toggle_homoquiralidade():
    global racemico
    # Alternar entre Homoquiral e Racêmico
    racemico = not racemico
    # Atualiza o texto do botão
    if not racemico:
        button_homoquiralidade.config(text="Homoquiral")
    else:
        button_homoquiralidade.config(text="Racêmico")

def toggle_laterais():
    global laterais
    # Alternar entre Homoquiral e Racêmico
    laterais = not laterais
    # Atualiza o texto do botão
    if not laterais:
        button_laterais.config(text="Sem Lig. Laterais")
    else:
        button_laterais.config(text="Com Lig. Laterais")

def toggle_hidropatica():
    global hidropatica
    # Alternar entre Homoquiral e Racêmico
    hidropatica = not hidropatica
    # Atualiza o texto do botão
    if hidropatica:
        button_hidropatica.config(text="Seq. Hidropáticas")
    else:
        button_hidropatica.config(text="Seq.  Aleatórias ")

def create_scrollable_frame(parent, title, data, entries, max_height=300):
    # Frame para a seção
    section_frame = tk.Frame(parent)
    section_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)

    # Título da seção
    label_header = tk.Label(section_frame, text=title, font=("Helvetica", 12, "bold"))
    label_header.pack(anchor="w")

    # Canvas para a barra de rolagem com altura máxima definida
    canvas = tk.Canvas(section_frame, width=200, height=max_height)
    scroll_y = tk.Scrollbar(section_frame, orient="vertical", command=canvas.yview)
    scroll_y.pack(side=tk.RIGHT, fill=tk.Y)

    # Frame interno do canvas
    input_frame = tk.Frame(canvas)
    canvas.create_window((0, 0), window=input_frame, anchor="nw")
    canvas.config(yscrollcommand=scroll_y.set)

    # Frame para os botões
    frame = tk.Frame(input_frame)
    frame.pack()

    # Botões de reset, zerar e distribuir
    reset_button = tk.Button(frame, text="Reset", command=lambda: reset_entries(entries, data), width=5)
    reset_button.pack(side=tk.RIGHT)

    zero_button = tk.Button(frame, text="Zerar", command=lambda: zero_entries(entries), width=5)
    zero_button.pack(side=tk.RIGHT)

    distr_button = tk.Button(frame, text="Distr.", command=lambda: distr_entries(entries), width=5)
    distr_button.pack(side=tk.RIGHT)

    # Adiciona os dados (aminoácidos)
    for amino_acido, values in data.items():
        full_name, proportion = values[:2]
        # Cria um frame para agrupar o label e o entry
        frame_entry = tk.Frame(input_frame)
        frame_entry.pack(anchor="w")

        # Adiciona o entry ao frame
        entry = tk.Entry(frame_entry, width=5)  # Ajuste a largura conforme necessário
        entry.insert(0, str(proportion))
        entry.pack(side=tk.LEFT)  # Adiciona o entry à direita do label

        # Adiciona o label ao frame
        label = tk.Label(frame_entry, text=f"{amino_acido} ({full_name}):")
        label.pack(side=tk.LEFT)  # Adiciona o label à esquerda do frame

        entries[amino_acido] = entry

    # Atualiza o tamanho do frame e do canvas
    input_frame.update_idletasks()
    canvas.config(scrollregion=canvas.bbox("all"))

    # Adiciona o canvas à seção
    canvas.pack(fill=tk.BOTH, expand=True)

    return entries


def reset_entries(entries, data):
    for amino_acido, entry in entries.items():
        entry.delete(0, tk.END)
        entry.insert(0, str(data[amino_acido][1]))  # Acessa o segundo elemento da tupla

def zero_entries(entries):
    for entry in entries.values():
        entry.delete(0, tk.END)
        entry.insert(0, '0')

def distr_entries(entries):
    for entry in entries.values():
        entry.delete(0, tk.END)
        entry.insert(0, 1)

# Inicializa a janela principal
root = tk.Tk()
root.title("Simulação de Polimerização by Lutzer")

# Define a janela sempre em primeiro plano
root.attributes('-topmost', True)

# Frame para o prompt
frame = tk.Frame(root)
frame.pack(side=tk.LEFT)

# Criando as seções com rolagem
entries_prot = create_scrollable_frame(frame, "AA's Padrão", amino_prot_proportion, {}, len(amino_prot_proportion)*10)
entries_nao_prot = create_scrollable_frame(frame, "AA's Alternativos", amino_nao_prot_proportion, {}, len(amino_nao_prot_proportion)*12)
entries_amines = create_scrollable_frame(frame, "Aminas Livres", amines_proportion, {}, len(amines_proportion)*20)

# Entradas para num_proteinas e tamanho_maximo - fora da rolagem
label_num_proteinas = tk.Label(root, text="Número de Proteínas:")
label_num_proteinas.pack(anchor="w", padx=10)
entry_num_proteinas = tk.Entry(root)
entry_num_proteinas.insert(0, "10000")
entry_num_proteinas.pack(padx=10)

label_tamanho_maximo = tk.Label(root, text="Tamanho Máximo:")
label_tamanho_maximo.pack(anchor="w", padx=10)
entry_tamanho_maximo = tk.Entry(root)
entry_tamanho_maximo.insert(0, "10")
entry_tamanho_maximo.pack(padx=10)

# Frame para o prompt
frame_opcoes = tk.Frame(root)
frame_opcoes.pack()

# Inicializa a variável homoquiral
racemico = True  # Começa como Homoquiral
# Botão que alterna entre Homoquiral e Racêmico
button_homoquiralidade = tk.Button(frame_opcoes, text="Racêmico", command=toggle_homoquiralidade)
button_homoquiralidade.pack(side=tk.LEFT, padx=10)

# Inicializa a variável ligações laterais
laterais = True  # Começa como Homoquiral
# Botão que alterna entre Homoquiral e Racêmico
button_laterais = tk.Button(frame_opcoes, text="Com Lig. Laterais", command=toggle_laterais)
button_laterais.pack(side=tk.LEFT, padx=10)

# Inicializa a variável ligações laterais
hidropatica = False  # Começa como Homoquiral
# Botão que alterna entre Homoquiral e Racêmico
button_hidropatica = tk.Button(frame_opcoes, text="Seq. Aleatórias", command=toggle_hidropatica)
button_hidropatica.pack(side=tk.LEFT, padx=10)

# Botão para calcular
calcular_button = tk.Button(root, text="Calcular Proporções", command=calcular_proporcoes)
calcular_button.pack(pady=10)

# Caixa de texto para resultados
result_box = scrolledtext.ScrolledText(root, width=50, height=25)

# Configurar a tag de negrito
result_box.tag_configure("normal", font=("Arial", 10))
result_box.tag_configure("title", font=("Arial", 12, "bold"))
result_box.tag_configure("bold", font=("Arial", 10, "bold"))
result_box.tag_configure("blue", foreground="blue", font=("Arial", 10))
result_box.tag_configure("blue under", foreground="blue", font=("Arial", 10, "underline"))
result_box.tag_configure("under", font=("Arial", 10, "underline"))
result_box.tag_configure("red", foreground="red", font=("Arial", 10, "bold"))
result_box.tag_configure('superscript', font=("Arial", 8), offset=5)

result_box.pack(padx=10, pady=10)

# Inicia a aplicação
root.mainloop()