# Distribution of number of tracks

The number of tracks is chosen randomly following this distribution

Tracks Number |  Probability
:------------:|:---------------:
   1          |     90   %
   2          |      9   %
   3          |      0.9 %
   4          |      0.1 %

but you can choose any number of tracks using the keyword `tracks_number` 
in `SimEvent`/`SimClunkyEvent` classes as follows
```python
event = SimEvent(tracks_number=7)
event = SimClunkyEvent(tracks_number=7)
```

This feature is implemented in `simulation/simulation.py`, `Simulate` class, 
parent of `SimEvent`/`SimClunkyEvent` classes.

