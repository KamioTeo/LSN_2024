from manim import *

# TODO: grafico con la dist che varia? aggiungo un punto alla volta

# Simulation folder
circle_path = "./data_pm20_pc80_I1k/Circle"
square_path = "./data_pm20_pc80_I1k/Square"

# Percorso posizioni cerchio
PATH_POS_C = circle_path + "/Positions.dat"
PATH_POS_S = square_path + "/Positions.dat"

# Percorso sequenze
PATH_SEQ_C = circle_path + "/minSequences.dat"
PATH_SEQ_S = square_path + "/minSequences.dat"


class DotPlacer(VMobject):
    """
    Mob che crea la forma esterna con la quale sono state disposte le città,
    e aggiunge un Dot per ognuna di esse. Ogni submobjects dall'1 in poi è un Dot
    in ordine di indice delle città, dalla prima alla 34
    """
    def __init__(self,
                 shape: VMobject,
                 dots_path: list,
                 scale_factor: float,
                 dots_color=BLUE,
                 **kwargs):
        super().__init__(**kwargs)
        self.shape = shape
        self.dots_path = dots_path
        self.scale_f = scale_factor
        self.dots_color = dots_color
        # creo la forma esterna
        self.add_body()
        # aggiungo i punti
        self.add_dots()

    def add_body(self):
        self.shape.scale(self.scale_f)
        self.add(self.shape)

    def add_dots(self):
        positions = self.positions_extractor(self.dots_path)
        dots = [Dot(self.scale_f*np.array(pos["coordinates"]), color=self.dots_color).set_z_index(2) for pos in positions]
        self.add(*dots)

    def positions_extractor(self, path: str) -> list:
        """
        Carico indice e coordinate per ogni città nelle due configurazioni, nel quadrato e sulla circonferenza.
        Li salvo in una lista di dizionari
        """
        with open(path) as posFile:
            rows = posFile.readlines()
            pos = []
            for row in rows[1:]:
                index, x, y = row.split()
                pos.append({"index": int(index), "coordinates": (float(x), float(y), 0)})
        return pos


class PathDrawer(VMobject):
    """
    Crea la spezzata che identifica una certa sequenza di città visitate
    (da dare in ingresso insieme ai dot, per sfruttare il loro posizionamento)
    """
    def __init__(self,
                 sequence: list,
                 shape: DotPlacer,
                 **kwargs):
        super().__init__(**kwargs)
        self.sequence = sequence
        # salvo i dots
        self.dots = shape[1:]
        self.add_path()

    def add_path(self):
        """Funzione che riordina la sequenza di dot come la sequenza passata"""
        new_dots = [self.dots[i-1] for i in self.sequence]
        self.set_points_as_corners(
            [d.get_center()
            for d in new_dots]
        )


def sequences_extractor(path: str) -> list:
    """
    Legge il file con le sequenze minime per generazione e restituisce una lista di dizionari,
    in cui in ognuno vi è la distanza e la sequenza.
    """
    with open(path) as seqFile:
        rows = seqFile.readlines()
        # contiene le sequenze di ogni popolazione, per ogni generazione (lista di liste di sequenze)
        gen_seq = []    
        for row in rows:
            dist, seq = row.split()
            seq_str = seq.split("-")
            seq_list = [int(index) for index in seq_str]
            gen_seq.append({"distance": float(dist), "sequence": seq_list})
    return gen_seq


class CircleAnimation(Scene):
    N_SEQ = 150 # numero sequenze da rappresentare
    def construct(self):
        # creo la forma e i punti di ogni città
        circle = DotPlacer(shape=Circle(radius=0.5, color=RED, stroke_width=5), dots_path=PATH_POS_C, scale_factor=5.5)
        self.play(Create(circle, lag_ratio=0.3, run_time=3))
        self.wait(0.5)
        # carico la lista di dizionari contenenti le distanze numeriche e le sequenze
        sequences_list = sequences_extractor(PATH_SEQ_C)
        # lista delle distanze numeriche per ogni sequenza
        distances = [
            sequences_list[i]["distance"]
            for i in range(0,self.N_SEQ)
        ]

        # creo gli assi
        ax= Axes(
                x_range=[-1,174,25], #START - END - STEP
                y_range=[-0.1,16,1],
                x_length=7, #lunghezza assi in unità manim
                y_length=6,
                axis_config={"include_numbers": True},
            )

        # affianco il cerchio e gli assi
        circle.generate_target()
        VGroup(circle.target, ax).arrange(RIGHT, buff=1)

        # creo le etichette del grafico (dopo aver spostato gli assi)
        x_label = ax.get_x_axis_label(
            Tex("Generations").scale(0.7), edge=DOWN, direction=DOWN, buff=0.5
        )
        y_label = ax.get_y_axis_label(
            Tex("Distance").scale(0.7).rotate(90 * DEGREES),
            edge=LEFT, #dice che voglio metterlo sull'asse parallelamente,
            direction=LEFT, #lo metto alla sinistra dell'asse
            buff=0.3,
        )
        # creo il testo con il valore numerico della distanza
        dist_tex = Tex("Distance: ").next_to(ax, UP, buff=-0.3).shift(LEFT)
        dists_dn = [
            DecimalNumber(distances[i], num_decimal_places=4, edge_to_fix=LEFT).next_to(dist_tex, RIGHT, buff=0.3)
            for i in range(0,self.N_SEQ)
        ]

        self.play(
            AnimationGroup(
                MoveToTarget(circle),
                AnimationGroup(
                    Create(ax),
                    Write(x_label),
                    Write(y_label),
                    lag_ratio=0.5,
                ),
                lag_ratio=0.6
            ),
            run_time=3
        )
        self.wait(0.5)

        # creo la lista con tutte le spezzate
        sequences = [
            PathDrawer(sequence=sequences_list[i]["sequence"], shape=circle, color=YELLOW)
            for i in range(0,self.N_SEQ)
        ]
        # creo la prima sequenza e scrivo le scritte
        self.play(
            AnimationGroup(
                Create(sequences[0]),
                AnimationGroup(
                    Write(dist_tex),
                    Write(dists_dn[0]),
                    lag_ratio=0.9,
                ),
                lag_ratio=0.8,    
                run_time=3,
                rate_func=linear
            )
        )
        self.wait()
        
        graph_dots = [Dot(ax.c2p(i, seq["distance"]), color=RED).scale(0.8) for i,seq in enumerate(sequences_list)]
        
        self.play(
            Succession(
                *[Transform(sequences[0], sequences[i+1])
                 for i in range(self.N_SEQ-1)],
                lag_ratio = 1.1
            ),
            Succession(
                *[Transform(dists_dn[0], dists_dn[i+1], rate_func=rate_functions.ease_in_expo)
                 for i in range(self.N_SEQ-1)],
                lag_ratio = 1.1
            ),
            Succession(
                *[Create(graph_dots[i])
                 for i in range(self.N_SEQ-1)],
                lag_ratio = 1.1
            ),
            run_time=20
        )
        self.wait()



class SquareAnimation(Scene):
    N_SEQ = 150 # numero sequenze da rappresentare
    def construct(self):
        # creo la forma e i punti di ogni città
        scale=5.5 # fattore di scala della figura quadrata e delle spezzate
        square = DotPlacer(shape=Square(side_length=1, color=RED, stroke_width=5), dots_path=PATH_POS_S, scale_factor=scale)
        # nel caso del quadrato le coordinate sono tutte positive (il centro è in (0.5,0.5))
        # mentre qua il centro è in (0,0) quindi devo traslare i punti in basso di (-0.5,-0.5) 
        # e moltiplicare per il fattore di scala
        square[1:].shift(scale*0.5*(DOWN+LEFT))
        self.play(Create(square, lag_ratio=0.3, run_time=3))
        self.wait(0.5)
        # carico la lista di dizionari contenenti le distanze numeriche e le sequenze
        sequences_list = sequences_extractor(PATH_SEQ_S)
        # lista delle distanze numeriche per ogni sequenza
        distances = [
            sequences_list[i]["distance"]
            for i in range(0,self.N_SEQ)
        ]

        # creo gli assi
        ax= Axes(
                x_range=[-1,174,25], #START - END - STEP
                y_range=[-0.1,16,1],
                x_length=7, #lunghezza assi in unità manim
                y_length=6,
                axis_config={"include_numbers": True},
            )

        # affianco il cerchio e gli assi
        square.generate_target()
        VGroup(square.target, ax).arrange(RIGHT, buff=1)

        # creo le etichette del grafico (dopo aver spostato gli assi)
        x_label = ax.get_x_axis_label(
            Tex("Generations").scale(0.7), edge=DOWN, direction=DOWN, buff=0.5
        )
        y_label = ax.get_y_axis_label(
            Tex("Distance").scale(0.7).rotate(90 * DEGREES),
            edge=LEFT, #dice che voglio metterlo sull'asse parallelamente,
            direction=LEFT, #lo metto alla sinistra dell'asse
            buff=0.3,
        )
        # creo il testo con il valore numerico della distanza
        dist_tex = Tex("Distance: ").next_to(ax, UP, buff=-0.3).shift(LEFT)
        dists_dn = [
            DecimalNumber(distances[i], num_decimal_places=4, edge_to_fix=LEFT).next_to(dist_tex, RIGHT, buff=0.3)
            for i in range(0,self.N_SEQ)
        ]

        self.play(
            AnimationGroup(
                MoveToTarget(square),
                AnimationGroup(
                    Create(ax),
                    Write(x_label),
                    Write(y_label),
                    lag_ratio=0.5,
                ),
                lag_ratio=0.6
            ),
            run_time=3
        )
        self.wait(0.5)

        # creo la lista con tutte le spezzate
        sequences = [
            PathDrawer(sequence=sequences_list[i]["sequence"], shape=square, color=YELLOW)
            for i in range(0,self.N_SEQ)
        ]
        # creo la prima sequenza e scrivo le scritte
        self.play(
            AnimationGroup(
                Create(sequences[0]),
                AnimationGroup(
                    Write(dist_tex),
                    Write(dists_dn[0]),
                    lag_ratio=0.9,
                ),
                lag_ratio=0.8,    
                run_time=3,
                rate_func=linear
            )
        )
        self.wait()
        
        graph_dots = [Dot(ax.c2p(i, seq["distance"]), color=RED).scale(0.8) for i,seq in enumerate(sequences_list)]
        
        self.play(
            Succession(
                *[Transform(sequences[0], sequences[i+1])
                 for i in range(self.N_SEQ-1)],
                lag_ratio = 1.1
            ),
            Succession(
                *[Transform(dists_dn[0], dists_dn[i+1], rate_func=rate_functions.ease_in_expo)
                 for i in range(self.N_SEQ-1)],
                lag_ratio = 1.1
            ),
            Succession(
                *[Create(graph_dots[i])
                 for i in range(self.N_SEQ-1)],
                lag_ratio = 1.1
            ),
            run_time=20
        )
        self.wait()
