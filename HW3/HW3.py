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
    def accelerate(self, amount):
        self.speed += amount
    def brake(self, amount):
        self.speed -= amount
        if self.speed < 0:
            self.speed = 0
    @classmethod
    def get_speed_limit(cls):
        return 60
    @staticmethod
    def honk():
        print("Honk honk!")
    @staticmethod
    def get_weather():
        response = requests.get('https://api.openmeteo.org')
        if response.status_code == 200:
            data = response.json()
            return data['weather']
        else:
            return "Unknown"
    @classmethod
    def get_num_cars_on_road(cls):
        return cls.num_cars_on_road
    @staticmethod
    def go_to_parking():
        Car.num_cars_on_road -= 1
    def __str__(self):
        return f"Car(speed={self.speed}, on_road={self.on_road})"
