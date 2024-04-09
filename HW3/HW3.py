import random
class IncreaseSpeedIterator:
    def __init__(self, car_instance, max_speed):
        self.car_instance = car_instance
        self.max_speed = max_speed
    def __iter__(self):
        return self
    def __next__(self):
        if self.car_instance.speed < self.max_speed:
            self.car_instance.speed += 1
            return self.car_instance.speed
        else:
            raise StopIteration
class DecreaseSpeedIterator:
    def __init__(self, car_instance):
        self.car_instance = car_instance
    def __iter__(self):
        return self
    def __next__(self):
        if self.car_instance.speed > 0:
            self.car_instance.speed -= 1
            return self.car_instance.speed
        else:
            raise StopIteration
class Car:
    num_cars_on_road = 0
    def __init__(self, speed=0, on_road=False):
        self.speed = speed
        self.on_road = on_road
        if self.on_road:
            Car.num_cars_on_road += 1
    def increase_speed(self, amount):
        self.speed += amount
    def decrease_speed(self, amount):
        self.speed -= amount
        if self.speed < 0:
            self.speed = 0
    def go_to_parking(self):
        self.on_road = False
        Car.num_cars_on_road -= 1
    def show_weather(self):
        weather_conditions = ['Sunny', 'Rainy', 'Cloudy', 'Stormy']
        return random.choice(weather_conditions)
# Using staticmethod
    def get_num_cars_on_road():
        return Car.num_cars_on_road
    def increase_speed_iterator(self, max_speed):
        return IncreaseSpeedIterator(self, max_speed)
    def decrease_speed_iterator(self):
        return DecreaseSpeedIterator(self)
