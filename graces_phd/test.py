import enum
from numbers import Number


class SnakeBreed(enum.Enum):
    SandBoa = "Sand Boa"
    BallPython = "Ball Python"


class Snake:
    def __init__(self, name: str, type: SnakeBreed) -> None:
        self.name = name
        self.type = type.value

    def sayHi(self):
        print("I'm " + self.name + " hisss")


def main():
    print("test")

    sum = addDaNums(1, 2)
    print(sum)

    riggs = Snake("Rigatoni", SnakeBreed.BallPython)
    ziti = Snake("Baby Snake", SnakeBreed.SandBoa)
    gracesSnakes: list[Snake] = [riggs, ziti]
    showMeSnakes(gracesSnakes)


def addDaNums(x: Number, y: Number):
    return (x+y)


def showMeSnakes(snakes: "list[Snake]"):
    for snake in snakes:
        print(snake.name)
        print(snake.type)
        snake.sayHi()
    return


main()
