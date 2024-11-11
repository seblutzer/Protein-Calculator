import random

import matplotlib.pyplot as plt



def wright_fisher_simulation(N, s, generations):

    """Simula a fixação de uma mutação benéfica em uma população usando o modelo de Wright-Fisher.

    Args:

        N (int): Tamanho da população.

        s (float): Coeficiente de seleção da mutação benéfica.

        generations (int): Número de gerações para simular.

    Returns:

        freq_history (list): Histórico da frequência da mutação ao longo das gerações.

    """

    freq = 1 / (2 * N)  # Frequência inicial da mutação em um indivíduo diploide

    freq_history = [freq]



    for _ in range(generations):

        # Probabilidade de seleção ajustada

        p = freq * (1 + s) / (freq * (1 + s) + (1 - freq))

        # Nova frequência após reprodução aleatória

        freq = random.binomial(2 * N, p) / (2 * N)

        freq_history.append(freq)

        # Verifica fixação ou perda

        if freq == 0 or freq == 1:

            break



    return freq_history



# Parâmetros da simulação

N = 1000          # Tamanho da população

s = 0.01          # Coeficiente de seleção

generations = 1000  # Número máximo de gerações



# Executa a simulação

freq_history = wright_fisher_simulation(N, s, generations)



# Plota os resultados

plt.plot(freq_history)

plt.xlabel('Gerações')

plt.ylabel('Frequência da Mutação')

plt.title('Simulação de Fixação de Mutação Benéfica')

plt.show()