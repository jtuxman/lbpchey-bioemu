# Guion — BioEmu y LBPCheY (8 min)

Ritmo de referencia: ~135 palabras/min (ritmo pausado de exposición). Texto hablado ≈ 1030 palabras ≈ 7:35; con pausas y transiciones ≈ 7:55.

| # | Diapositiva | Duración | Acumulado |
|---|---|---|---|
| 1 | Portada | 0:20 | 0:20 |
| 2 | La función depende de un ensamble | 0:50 | 1:10 |
| 3 | Evoformer + difusión, en tres etapas | 0:55 | 2:05 |
| 4 | Recupera estados funcionales | 0:50 | 2:55 |
| 5 | Reproduce paisajes de MD | 0:40 | 3:35 |
| 6 | PPFT | 0:45 | 4:20 |
| 7 | LBPCheY: dos folds y una bisagra | 1:00 | 5:20 |
| 8 | Un estado dominante, centrado en AF3 | 1:00 | 6:20 |
| 9 | ¿Qué podemos afirmar con 246 muestras? | 0:50 | 7:10 |
| 10 | Conclusiones | 0:45 | 7:55 |

> Si vas tarde al llegar a la diapositiva 7 (más de 4:35), salta los textos marcados **[opcional]**.

---

## 1 · Portada — 0:20 (0:00 → 0:20)

Hoy les presento BioEmu, un modelo generativo publicado en Science este año que predice no una estructura, sino el ensamble de equilibrio de una proteína. Primero veremos qué hace y qué tan bien funciona, y después lo aplicaremos a nuestra quimera, LBPCheY.

---

## 2 · La función depende de un ensamble — 0:50 (0:20 → 1:10)

*(Señalar panel A)* Las proteínas no son estáticas. La actina, por ejemplo, funciona porque cicla entre estados abierto y cerrado, y cada uno de esos pasos es un equilibrio entre conformaciones.

*(Señalar las tres tarjetas de abajo)* Hoy tenemos dos herramientas para esto, y ninguna es suficiente. AlphaFold nos da una sola estructura: una fotografía, sin poblaciones. La dinámica molecular sí nos da el ensamble, pero para cruzar barreras de energía necesita microsegundos o milisegundos de simulación, que son meses o años de GPU.

*(Señalar panel B)* BioEmu propone un atajo: a partir de la secuencia, generar muestras directamente de la distribución de equilibrio. Miles de estructuras por hora en una sola GPU.

---

## 3 · Evoformer + difusión, en tres etapas — 0:55 (1:10 → 2:05)

*(Señalar Fig. 1C)* La arquitectura tiene dos partes. Primero, el Evoformer de AlphaFold2 codifica la secuencia y su alineamiento múltiple. Esa información condiciona a un modelo de difusión, que parte de ruido y, paso a paso, lo convierte en una estructura. Cada vez que corremos el modelo obtenemos una conformación distinta.

*(Señalar Fig. 1D)* Lo realmente nuevo es cómo lo entrenan, en tres etapas. Primero, con unos doscientos millones de estructuras de la base de AlphaFold, para aprender diversidad estructural. Después, con más de doscientos milisegundos de dinámica molecular, para aprender la física del equilibrio. Y por último, con medio millón de mediciones experimentales de estabilidad, para calibrar las poblaciones.

La idea clave es que integra tres tipos de datos muy distintos: estructuras, simulaciones y mediciones experimentales.

---

## 4 · Recupera estados funcionales — 0:50 (2:05 → 2:55)

La primera pregunta es si BioEmu encuentra conformaciones que sabemos que existen. *(Señalar los mapas)* Aquí vemos tres proteínas con movimientos grandes de dominio. En verde están las muestras de BioEmu y en gris las estructuras del PDB. Los mapas muestran que el modelo puebla tanto el estado abierto como el cerrado.

*(Señalar LAO-binding protein)* Fíjense en esta: es una proteína periplásmica de unión, de la misma familia que LBP. Nos dice que BioEmu ya captura el abrir y cerrar de este tipo de proteínas.

*(Señalar los números)* En total recupera el 83 % de los movimientos de dominio y el 86 % de los bolsillos crípticos. Además, el resultado no depende de qué tan parecida sea la proteína al set de entrenamiento.

---

## 5 · Reproduce paisajes de MD — 0:40 (2:55 → 3:35)

La segunda pregunta es si, además de encontrar los estados, les asigna las poblaciones correctas. *(Señalar los mapas MD vs BioEmu)* Para las proteínas de plegamiento rápido de DESRES, el paisaje de BioEmu reproduce los mismos mínimos que simulaciones de milisegundos, con un error de 0.74 kcal/mol, del orden del error experimental.

*(Señalar la gráfica de costo)* Y lo hace con unas diez mil veces menos horas de GPU.

**[opcional]** Lo mismo se cumple en dominios CATH, en una proteína desordenada y en una tetraspanina de membrana.

---

## 6 · PPFT — 0:45 (3:35 → 4:20)

¿Cómo usan los datos experimentales si no hay estructuras de cada mutante? Con una técnica llamada PPFT. *(Recorrer el esquema de izquierda a derecha)* Generan muestras rápidas, clasifican cada una como plegada o desplegada, y la proporción entre ambas da una energía libre. Esa energía se compara con el valor experimental, y el error se usa para corregir el modelo.

*(Señalar las gráficas y números)* El resultado: menos de 0.9 kcal/mol de error en estabilidad y menos de 0.8 en el efecto de mutaciones.

Y esto es justo lo que usaremos ahora: en BioEmu, contar cuántas estructuras caen en un estado equivale a medir su energía libre.

---

## 7 · LBPCheY: dos folds y una bisagra — 1:00 (4:20 → 5:20)

Pasamos a nuestro experimento. LBPCheY es una quimera de 339 residuos: el dominio CheY está insertado en medio de la proteína de unión a leucina, LBP. Sabemos que los dos dominios se pliegan bien. La pregunta es si se mueven uno respecto al otro.

*(Señalar el pipeline)* Partimos del modelo de AlphaFold3, con un pTM de 0.81, generamos 300 muestras con BioEmu, quedaron 246 válidas, y les reconstruimos las cadenas laterales.

*(Señalar los histogramas)* Al medir el RMSD por dominio, CheY se mantiene muy cerca del modelo, en 0.8 angstroms: la inserción no lo desestabiliza. LBP es más variable, en parte porque está partida en dos segmentos y tiene su propia bisagra.

Lo importante es que el RMSD global es mayor que el de cada dominio. Eso significa que los dominios no se están desplegando: se están moviendo uno respecto al otro.

---

## 8 · Un estado dominante, centrado en AF3 — 1:00 (5:20 → 6:20)

*(Señalar panel b, derecha)* Para ver ese movimiento construimos la superficie de energía libre contando muestras. En el eje horizontal está la distancia entre los centros de los dos dominios, y en el vertical cuánto gira CheY respecto al modelo de AlphaFold3.

Vemos tres cosas. Primero, hay un solo pozo principal, y la estrella, que es AlphaFold3, cae justo en él. La predicción estática corresponde al estado más probable.

Segundo, ese pozo es ancho: CheY puede girar unos cincuenta grados hacia cada lado y la distancia varía unos ocho angstroms sin cruzar barreras. No son dos estados discretos: es una bisagra flexible.

Tercero, alrededor aparecen islas a unas 3 kcal/mol, con giros mucho mayores o con los dominios más separados. Son accesibles térmicamente, pero cada una está formada por solo una o dos muestras.

---

## 9 · ¿Qué podemos afirmar con 246 muestras? — 0:50 (6:20 → 7:10)

*(Señalar perfil izquierdo)* En los perfiles de una dimensión, la distancia muestra un pozo plano y asimétrico: separar los dominios cuesta menos que acercarlos. *(Señalar perfil derecho)* El giro tiene su mínimo exactamente en cero, que es AlphaFold3.

*(Señalar la tarjeta verde)* Lo que sí podemos afirmar: hay una bisagra flexible alrededor del modelo de AlphaFold3, y las 246 muestras tienen ambos dominios plegados.

*(Señalar la tarjeta roja)* Y lo que debemos tomar con cautela: con 246 muestras, la energía libre máxima que podemos resolver es de unas 3.3 kcal/mol. Por eso los pozos secundarios, justo en ese límite, son hipótesis y no resultados. Además, BioEmu tiene un error de alrededor de 1 kcal/mol, no incluye ligandos y no nos da cinética.

---

## 10 · Conclusiones — 0:45 (7:10 → 7:55)

En resumen: BioEmu genera ensambles de equilibrio a partir de la secuencia, recupera estados funcionales y reproduce paisajes de dinámica molecular con un error de alrededor de 1 kcal/mol, a una fracción del costo. No reemplaza a la dinámica molecular, la complementa.

En LBPCheY nos mostró, en horas, que ambos dominios se conservan y que están unidos por una bisagra flexible con un estado dominante que coincide con AlphaFold3.

Los siguientes pasos son generar más muestras para confirmar los estados raros y usar dinámica molecular para ver qué hace la leucina con esta bisagra. Gracias.

---

### Respuestas cortas para preguntas probables

- **¿Por qué solo 246 de 300?** BioEmu descarta muestras con choques estéricos o geometría inválida del backbone.
- **¿Qué es TIC en la Fig. 3?** Componentes independientes en el tiempo: las coordenadas más lentas de la simulación.
- **¿Cómo calculas ΔG?** ΔG = −kT·ln(pᵢ/p_max), con T = 300 K (kT = 0.596 kcal/mol), contando muestras por bin.
- **¿Por qué no modelaron la leucina con BioEmu?** BioEmu solo modela una cadena proteica, sin ligandos.
- **¿Cuántas muestras harían falta?** Para resolver ~4 kcal/mol, del orden de 1000 o más (kT·ln(1000) ≈ 4.1 kcal/mol).
