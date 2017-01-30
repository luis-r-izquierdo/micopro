;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; GNU GENERAL PUBLIC LICENSE ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; MiCoPro (Mixing in Contagion Processes) is a model designed to analyze
;; the consequences of mixing two groups that have different predispositions
;; to adopt a certain trait (or infection).
;; Copyright (C) 2017 Luis R. Izquierdo, Segismundo S. Izquierdo & Dunia López-Pintado
;;
;; This program is free software: you can redistribute it and/or modify
;; it under the terms of the GNU General Public License as published by
;; the Free Software Foundation, either version 3 of the License, or
;; (at your option) any later version.
;;
;; This program is distributed in the hope that it will be useful,
;; but WITHOUT ANY WARRANTY; without even the implied warranty of
;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;; GNU General Public License for more details.
;;
;; You should have received a copy of the GNU General Public License
;; along with this program.  If not, see <http://www.gnu.org/licenses/>.
;;
;; Contact information:
;; Luis R. Izquierdo
;;   University of Burgos, Spain.
;;   e-mail: lrizquierdo@ubu.es


;;;;;;;;;;;;;;;;;
;;; Variables ;;;
;;;;;;;;;;;;;;;;;


globals [

  pi11 pi22
  probability-of-recovery
  exp-probability-of-infection-g1
  exp-probability-of-infection-g2
  avg-probability-of-infection-g1
  avg-probability-of-infection-g2
  fraction-of-interacting-agents

  avg-lambda-g1
  avg-lambda-g2

  ;; Variables for the agent-based model

  n-of-g1 ;; number of agents in group 1
  n-of-g2 ;; number of agents in group 2

  g1 ;; agentset containing group 1 agents
  g2 ;; agentset containing group 2 agents

  %-infected-agents-g1-ab
  %-infected-agents-g2-ab
  %-infected-agents-ab

  list-of-latest-%-infected-agents-ab
  list-of-latest-%-infected-agents-g1-ab
  list-of-latest-%-infected-agents-g2-ab

  ;; Variables for the mean dynamics

  freq-of-infected-agents-g1
  freq-of-infected-agents-g2

  %-infected-agents-g1-md
  %-infected-agents-g2-md
  %-infected-agents-md

  list-of-stable-%-infected-agents-md

  critical-mixing

  steady-state?
  %-infected-agents-x-ticks-ago %-infected-agents-now
]

breed [individuals individual]

individuals-own [
  my-group
  mate

  infected?
  infected?-next-period

  probability-of-infection
]

;;;;;;;;;;;;;;;;;;;;;;;;
;;; Setup Procedures ;;;
;;;;;;;;;;;;;;;;;;;;;;;;

to startup
  clear-all
  no-display
  setup-specific-variables
  setup-agents
  setup-variables
  do-initial-infection
  gather-data-ab
  gather-data-md
  reset-ticks
end

to setup-specific-variables
  set pi11 (1 - mixing / 100)
  set pi22 (1 - mixing / 100)
  set probability-of-recovery 0.03125
    ;; Assuming fraction-of-interacting-agents = 0.5, the upper line
    ;; imposes an upper limit on the value of lambda: lambda <= 16.
    ;; Since there could be variability, the population expected value of lambda
    ;; (the value that is set in the interface) should be <= 8.
    ;; I.e. if the user sets lamda = 8, there could an agent with %-variability = 100,
    ;; who would have a value of lambda = 16, and then probability-of-infection = 1.
  set fraction-of-interacting-agents 0.5
  set exp-probability-of-infection-g1 (expected-lambda-group-1 * probability-of-recovery / fraction-of-interacting-agents)
  set exp-probability-of-infection-g2 (expected-lambda-group-2 * probability-of-recovery / fraction-of-interacting-agents)
end

to setup-agents

  set n-of-g1 round (n-of-agents-in-total * ifelse-value (pi11 = pi22)
    [0.5] [ (1 - pi22) / (2 - pi11 - pi22)])
    ;; done to avoid division by zero when pi11 = pi22 = 1

  ;; group 1 agents
  create-individuals round (n-of-g1 / 2) [
    set my-group 1
    set mate nobody
    let random-% (-1 + random-float 2) * %-variability
    set probability-of-infection (exp-probability-of-infection-g1 * (1 + random-% / 100))
    set hidden? true
    hatch-individuals 1 [set probability-of-infection (exp-probability-of-infection-g1 * (1 - random-% / 100))]
  ]

  ;; group 2 agents
  set n-of-g2 (n-of-agents-in-total - n-of-g1)
  create-individuals round (n-of-g2 / 2) [
    set my-group 2
    set mate nobody
    let random-% (-1 + random-float 2) * %-variability
    set probability-of-infection (exp-probability-of-infection-g2 * (1 + random-% / 100))
    set hidden? true
    hatch-individuals 1 [set probability-of-infection (exp-probability-of-infection-g2 * (1 - random-% / 100))]
  ]

end

to setup-variables

  ;; Variables for the agent-based model

  set g1 individuals with [my-group = 1]
  set g2 individuals with [my-group = 2]

  set list-of-latest-%-infected-agents-ab []
  set list-of-latest-%-infected-agents-g1-ab []
  set list-of-latest-%-infected-agents-g2-ab []

  ;; Variables for the mean dynamics

  set freq-of-infected-agents-g1 initial-%-of-infected / 100
  set freq-of-infected-agents-g2 initial-%-of-infected / 100

  set list-of-stable-%-infected-agents-md []

  set avg-probability-of-infection-g1 mean [probability-of-infection] of g1
  set avg-probability-of-infection-g2 mean [probability-of-infection] of g2

  set avg-lambda-g1 avg-probability-of-infection-g1 * fraction-of-interacting-agents / probability-of-recovery
  set avg-lambda-g2 avg-probability-of-infection-g2 * fraction-of-interacting-agents / probability-of-recovery

  update-critical-mixing

  set steady-state? false

end

to update-critical-mixing
  set critical-mixing 100 * (1 - ifelse-value (avg-lambda-g1 + avg-lambda-g2 - 2 * avg-lambda-g1 * avg-lambda-g2 = 0)
    ["N/A"] [(1 - avg-lambda-g1 * avg-lambda-g2) / (avg-lambda-g1 + avg-lambda-g2 - 2 * avg-lambda-g1 * avg-lambda-g2)])
end

to do-initial-infection
  ask individuals [set infected? false]
  ask (turtle-set
    (n-of (round (n-of-g1 * initial-%-of-infected / 100)) g1)
    (n-of (round (n-of-g2 * initial-%-of-infected / 100)) g2)) [
      set infected? true
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Run-time procedures ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;

to go

  if agent-based-model? [
    if any? individuals with [infected?] [
      make-couples
      ask individuals [do-infections-and-cures]
      ask individuals [set infected? infected?-next-period]
    ]
    gather-data-ab
  ]

  ;; Mean dynamics

  repeat time-resolution [update-mean-dynamics]
  gather-data-md

  ;; Look for steady points in the mean dynamics
  if (ticks mod 50 = 0) [check-stability]

   ;; END Mean dynamics

  if ticks > 1000 [
    set list-of-latest-%-infected-agents-ab lput %-infected-agents-ab list-of-latest-%-infected-agents-ab
    set list-of-latest-%-infected-agents-g1-ab lput %-infected-agents-g1-ab list-of-latest-%-infected-agents-g1-ab
    set list-of-latest-%-infected-agents-g2-ab lput %-infected-agents-g2-ab list-of-latest-%-infected-agents-g2-ab
    set list-of-stable-%-infected-agents-md lput %-infected-agents-md list-of-stable-%-infected-agents-md
  ]

  tick

end

to make-couples
  ;; The following code is written here to allow the user
  ;; to change the value of mixing at runtime. For that, we may need
  ;; to adjust the number of agents.

  set pi11 (1 - mixing / 100)
  set pi22 (1 - mixing / 100)

  adjust-num-agents

  ;; We compute the following numbers here
  ;; (rather than just once at the beginning of the simulation)
  ;; to allow the user modify the value of the fraction-of-interacting-agents.

  ;; the following agentsets contain an even number of agents
  let g1-to-interact-with-g1 n-of (2 * floor (fraction-of-interacting-agents * n-of-g1 * pi11 / 2)) g1
  let g2-to-interact-with-g2 n-of (2 * floor (fraction-of-interacting-agents * n-of-g2 * pi22 / 2)) g2

  ask individuals [set mate nobody]

  pair-yourselves g1-to-interact-with-g1
  pair-yourselves g2-to-interact-with-g2

  let n-of-mixed-interactions floor (fraction-of-interacting-agents * n-of-g1 * (1 - pi11))
    ;; or floor (fraction-of-interacting-agents * n-of-g2 * (1 - pi22))
  let g1-to-interact-with-g2 n-of n-of-mixed-interactions (g1 with [mate = nobody])
  let g2-to-interact-with-g1 n-of n-of-mixed-interactions (g2 with [mate = nobody])

  let g1-to-interact-with-g2-list shuffle sort g1-to-interact-with-g2
  let g2-to-interact-with-g1-list sort g2-to-interact-with-g1
  (foreach g1-to-interact-with-g2-list g2-to-interact-with-g1-list [ [?1 ?2] ->
    ask ?1 [set mate ?2]
    ask ?2 [set mate ?1]
  ])

  ;; the following alternative code is less efficient but easier to understand:
  ;  ask g1-to-interact-with-g2 [
  ;    set mate one-of (g2-to-interact-with-g1 with [mate = nobody])
  ;    ask mate [set mate myself]
  ;  ]

  ;; normally, there can be up to 2 agents of each group left out even if fraction-of-interacting-agents = 1.
  ;; (fraction-of-interacting-agents * n-of-g1 * pi11) could be something like 3.8, which would lead to
  ;; g1-to-interact-with-g1 = 2

end

to adjust-num-agents

  ;; group 1 agents
  set n-of-g1 round (n-of-agents-in-total * ifelse-value (pi11 = pi22)
    [0.5] [ (1 - pi22) / (2 - pi11 - pi22)])
    ;; done to avoid division by zero when pi11 = pi22 = 1

  let adjustment-g1 (n-of-g1 - (count g1))

  if adjustment-g1 != 0 [
    ifelse adjustment-g1 > 0
    [
      set exp-probability-of-infection-g1 (expected-lambda-group-1 * probability-of-recovery / fraction-of-interacting-agents)

      create-individuals adjustment-g1 [
        set my-group 1
        set mate nobody
        set probability-of-infection (exp-probability-of-infection-g1 * (1 + (-1 + random-float 2) * %-variability / 100))
        set hidden? true
        set infected? ifelse-value (random-float 100.0 < %-infected-agents-g1-ab) [true][false]
      ]
    ]
    [
      ask n-of (0 - adjustment-g1) g1 [
        if (mate != nobody) [ ask mate [set mate nobody] ]
        die
      ]
    ]

    set g1 individuals with [my-group = 1]
    set avg-probability-of-infection-g1 mean [probability-of-infection] of g1
    set avg-lambda-g1 avg-probability-of-infection-g1 / probability-of-recovery
  ]

  ;; group 2 agents

  set n-of-g2 (n-of-agents-in-total - n-of-g1)
  let adjustment-g2 (n-of-g2 - (count g2))

  if adjustment-g2 != 0 [
    ifelse adjustment-g2 > 0
    [
      set exp-probability-of-infection-g2 (expected-lambda-group-2 * probability-of-recovery / fraction-of-interacting-agents)
      create-individuals adjustment-g2 [
        set my-group 2
        set mate nobody
        set probability-of-infection (exp-probability-of-infection-g2 * (1 + (-1 + random-float 2) * %-variability / 100))
        set hidden? true
        set infected? ifelse-value (random-float 100.0 < %-infected-agents-g2-ab) [true][false]
      ]
    ]
    [
      ask n-of (0 - adjustment-g2) g2 [
        if (mate != nobody)  [ ask mate [set mate nobody] ]
        die
      ]
    ]

    set g2 individuals with [my-group = 2]
    set avg-probability-of-infection-g2 mean [probability-of-infection] of g2
    set avg-lambda-g2 avg-probability-of-infection-g2 / probability-of-recovery
  ]

  update-critical-mixing
end

to pair-yourselves [agent-set]
  let n count agent-set
  let agent-list shuffle sort agent-set
  let first-half sublist agent-list 0 (n / 2)
  let second-half sublist agent-list (n / 2) n
  (foreach first-half second-half [ [?1 ?2] ->
    ask ?1 [set mate ?2]
    ask ?2 [set mate ?1]
  ])

  ;; the following alternative code is less efficient but easier to understand:
  ;  ask agent-set [
  ;    if (mate = nobody) [
  ;      set mate one-of (agent-set with [mate = nobody and self != myself])
  ;      ask mate [set mate myself]
  ;    ]
  ;  ]
end

to do-infections-and-cures
  ;; this procedure only changes the value of infected?-next-period
  ;; it does not change the value of infected?
  ;; this corresponds to synchronous updating

  set infected?-next-period infected?

  ;; note that you only enter one of the branches of the code;
  ;; which one depends on whether you are infected or not
  ifelse infected?
  [
    if random-float 1.0 < probability-of-recovery [set infected?-next-period false]
  ]
  [
    if (mate != nobody) [
      if random-float 1.0 < probability-of-infection [
        if ([infected?] of mate) [set infected?-next-period true]
      ]
    ]
  ]
end

to update-mean-dynamics

  set pi11 (1 - mixing / 100)
  set pi22 (1 - mixing / 100)

  let prob-meeting-infected-agent-g1 fraction-of-interacting-agents *
                                     (pi11 * freq-of-infected-agents-g1 + (1 - pi11) * freq-of-infected-agents-g2)
  let prob-meeting-infected-agent-g2 fraction-of-interacting-agents *
                                     ((1 - pi22) * freq-of-infected-agents-g1 + pi22 * freq-of-infected-agents-g2)

  ;; rate at which a susceptible agent of group i becomes infected
  let infection-rate-g1 avg-probability-of-infection-g1 * prob-meeting-infected-agent-g1
  let infection-rate-g2 avg-probability-of-infection-g2 * prob-meeting-infected-agent-g2

  ;; it is best to choose powers of 2 for the value of time-resolution (since we are dividing)
  set freq-of-infected-agents-g1 freq-of-infected-agents-g1 + (((1 - freq-of-infected-agents-g1) * infection-rate-g1 - freq-of-infected-agents-g1 * probability-of-recovery) / time-resolution)
  set freq-of-infected-agents-g2 freq-of-infected-agents-g2 + (((1 - freq-of-infected-agents-g2) * infection-rate-g2 - freq-of-infected-agents-g2 * probability-of-recovery) / time-resolution)

end

;;;;;;;;;;;;;;;;;;;;;;;;
;;;    Statistics    ;;;
;;;;;;;;;;;;;;;;;;;;;;;;

to gather-data-ab

  ;; Agent-based model
  set %-infected-agents-g1-ab 100 * (count g1 with [infected?]) / n-of-g1
  set %-infected-agents-g2-ab 100 * (count g2 with [infected?]) / n-of-g2
  set %-infected-agents-ab (n-of-g1 * %-infected-agents-g1-ab + n-of-g2 * %-infected-agents-g2-ab) / (n-of-g1 + n-of-g2)

end

to gather-data-md

  ;; Mean dynamics
  set %-infected-agents-g1-md 100 * freq-of-infected-agents-g1
  set %-infected-agents-g2-md 100 * freq-of-infected-agents-g2
  set %-infected-agents-md (%-infected-agents-g1-md * n-of-g1 + %-infected-agents-g2-md * n-of-g2) / (n-of-g1 + n-of-g2)

end

to check-stability
  let %-infected-agents-2x-ticks-ago %-infected-agents-x-ticks-ago
  set %-infected-agents-x-ticks-ago %-infected-agents-now
  set %-infected-agents-now %-infected-agents-md

  set steady-state? not-significantly-different?
    (list %-infected-agents-now %-infected-agents-x-ticks-ago %-infected-agents-2x-ticks-ago)
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Supporting procedures ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report not-significantly-different? [l]
  while [length l > 1] [
    let first-item first l
    let rest but-first l
    let all-similar? reduce [ [?1 ?2] -> ?1 and ?2 ] (map [ [?1] -> abs (?1 - first-item) < 0.001 ] rest)
    if not all-similar? [report false]
    set l rest
  ]
  report true
end
@#$#@#$#@
GRAPHICS-WINDOW
286
12
459
186
-1
-1
82.5
1
10
1
1
1
0
1
1
1
0
1
0
1
0
0
1
ticks
30.0

SLIDER
27
313
278
346
initial-%-of-infected
initial-%-of-infected
0
100
1.0
0.1
1
%
HORIZONTAL

BUTTON
24
10
90
43
setup
startup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
100
10
181
43
go once
go
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
0

BUTTON
190
10
278
43
NIL
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
0

PLOT
286
10
570
223
% Infected - Agent-based model
NIL
NIL
0.0
10.0
0.0
0.0
true
true
"" ""
PENS
"g1" 1.0 0 -13345367 true "" "plot %-infected-agents-g1-ab"
"g2" 1.0 0 -2674135 true "" "plot %-infected-agents-g2-ab"
"all" 1.0 0 -7500403 true "" "plot %-infected-agents-ab"

PLOT
287
228
569
433
% infected - Mean dynamics
NIL
NIL
0.0
10.0
0.0
0.0
true
true
"" ""
PENS
"g1" 1.0 0 -13345367 true "" "plot %-infected-agents-g1-md"
"g2" 1.0 0 -2674135 true "" "plot %-infected-agents-g2-md"
"all" 1.0 0 -7500403 true "" "plot %-infected-agents-md"

SLIDER
26
399
280
432
time-resolution
time-resolution
1
128
128.0
1
1
NIL
HORIZONTAL

MONITOR
27
350
159
395
critical mixing (%)
critical-mixing
2
1
11

MONITOR
577
46
711
91
%-infected-g1-abm
%-infected-agents-g1-ab
2
1
11

MONITOR
577
94
711
139
%-infected-g2-abm
%-infected-agents-g2-ab
2
1
11

MONITOR
575
229
710
274
%-infected-g1-md
%-infected-agents-g1-md
2
1
11

MONITOR
575
277
710
322
%-infected-g2-md
%-infected-agents-g2-md
2
1
11

MONITOR
574
388
670
433
NIL
steady-state?
0
1
11

MONITOR
204
52
279
97
time
ticks
17
1
11

SWITCH
577
10
763
43
agent-based-model?
agent-based-model?
0
1
-1000

MONITOR
575
325
734
370
NIL
%-infected-agents-md
2
1
11

MONITOR
577
142
736
187
%-infected-agents-abm
%-infected-agents-ab
2
1
11

SLIDER
24
57
200
90
mixing
mixing
0
100
50.0
1
1
%
HORIZONTAL

SLIDER
25
104
278
137
n-of-agents-in-total
n-of-agents-in-total
1000
10000
5000.0
1000
1
NIL
HORIZONTAL

SLIDER
25
148
277
181
expected-lambda-group-1
expected-lambda-group-1
0
6
0.5
0.05
1
NIL
HORIZONTAL

SLIDER
25
185
277
218
expected-lambda-group-2
expected-lambda-group-2
0
6
2.0
0.05
1
NIL
HORIZONTAL

SLIDER
25
228
278
261
%-variability
%-variability
0
100
0.0
1
1
%
HORIZONTAL

MONITOR
26
265
150
310
avg lambda-g1
avg-lambda-g1
3
1
11

MONITOR
159
265
278
310
avg lambda-g2
avg-lambda-g2
3
1
11

@#$#@#$#@
## WHAT IS IT?

(a general understanding of what the model is trying to show or explain)

## HOW IT WORKS

(what rules the agents use to create the overall behavior of the model)

## HOW TO USE IT

(how to use the model, including a description of each of the items in the Interface tab)

## THINGS TO NOTICE

(suggested things for the user to notice while running the model)

## THINGS TO TRY

(suggested things for the user to try to do (move sliders, switches, etc.) with the model)

## EXTENDING THE MODEL

(suggested things to add or change in the Code tab to make the model more complicated, detailed, accurate, etc.)

## NETLOGO FEATURES

(interesting or unusual features of NetLogo that the model uses, particularly in the Code tab; or where workarounds were needed for missing features)

## RELATED MODELS

(models in the NetLogo Models Library and elsewhere which are of related interest)

## CREDITS AND REFERENCES

(a reference to the model's URL on the web if it has one, as well as any other necessary credits, citations, and links)
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.0
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="experiment" repetitions="10" runMetricsEveryStep="false">
    <setup>startup</setup>
    <go>go</go>
    <timeLimit steps="2001"/>
    <metric>(length list-of-latest-%-infected-agents-ab)</metric>
    <metric>(mean list-of-latest-%-infected-agents-t1-ab)</metric>
    <metric>%-infected-agents-t1-md</metric>
    <metric>(mean list-of-latest-%-infected-agents-t2-ab)</metric>
    <metric>%-infected-agents-t2-md</metric>
    <metric>(mean list-of-latest-%-infected-agents-ab)</metric>
    <metric>%-infected-agents-md</metric>
    <metric>steady-state?</metric>
    <enumeratedValueSet variable="n-of-agents-in-total">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="lambda-t1">
      <value value="0.5"/>
      <value value="2"/>
      <value value="2.5"/>
      <value value="4"/>
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="lambda-t2">
      <value value="0.25"/>
      <value value="0.5"/>
      <value value="1"/>
      <value value="2"/>
      <value value="2.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="time-resolution">
      <value value="128"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="agent-based-model?">
      <value value="true"/>
    </enumeratedValueSet>
    <steppedValueSet variable="π" first="0" step="0.05" last="1"/>
    <enumeratedValueSet variable="initial-%-of-infected">
      <value value="50"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
