from enum import Enum

class Translate(Enum):
    LEFT        = 2
    RIGHT       = 1 
    FRONT       = 5
    BACKWARDS   = 6
    UP          = 4
    DOWN        = 3
    
class Scale(Enum):
    ZOOMIN        = 8
    ZOOMOUT       = 9 
    
class Rotate(Enum): 
    # based on Right-Hand Coordinate 
    # the variable  behind '_' is the normal vector
    # e.g. "CLOCKWISE_Z" means clockwise rotation relative to Z axis
    
    #    → → →
    #    ↑   ↓
    #    ← ← ←
    CLOCKWISE_Z           = 10
    #    ← ← ←
    #    ↓   ↑
    #    → → →
    COUNTERCLOCKWISE_Z    = 11
    
    CLOCKWISE_X           = 20
    
    COUNTERCLOCKWISE_X    = 21
    
    CLOCKWISE_Y           = 30
    
    COUNTERCLOCKWISE_Y    = 31
    
    FLIP_HORIZONTAL       = 56
    
