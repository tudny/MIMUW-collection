import argparse
import json
import logging.handlers
import random
import threading
import time
import logging
from kafka import KafkaProducer
from typing import Final

LOG_FILE: Final[str] = "data_generator.log"

logging.basicConfig(
    filename=LOG_FILE,
    filemode="a",
    format="%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s",
    datefmt="%H:%M:%S",
    level=logging.WARN,
)

handler = logging.handlers.RotatingFileHandler(
    LOG_FILE,
    maxBytes=1000000,
    backupCount=5,
)

LOGGER = logging.getLogger(__name__)

STOP_EVENT = threading.Event()


def stream_1_generator():
    index = 0
    while True:
        index += 1
        if random.uniform(0, 1) < 0.05:
            continue
        yield (index, random.uniform(0, 1))


def stream_2_generator():
    break_point = 300_000
    index = 0
    while index < break_point:
        index += 1
        if random.uniform(0, 1) < 0.1:
            continue
        yield (index, random.gauss(2, 1))
    while index < break_point * 2:
        index += 1
        if random.uniform(0, 1) < (index - break_point) / break_point:
            yield (index, random.uniform(0, 1))
        else:
            yield (index, random.gauss(2, 1))

    while True:
        index += 1
        if random.uniform(0, 1) < 0.01:
            continue
        yield (index, random.uniform(0, 1))


def producer_job(topic: str, stream_generator, server: str):
    producer = KafkaProducer(
        bootstrap_servers=server,
        value_serializer=lambda v: json.dumps(v).encode("utf-8"),
    )
    gen = stream_generator()
    logging.info("Starting producer for %s", topic)
    logging.info("Server: %s", server)
    while not STOP_EVENT.is_set():
        time_point, value = next(gen)
        # print(f"Sending {time_point} {value} to {topic}")
        logging.info("Sending %s %s to %s", time_point, value, topic)
        producer.send(topic, {"time_point": time_point, "value": value})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--topic1", type=str, default="topic1")
    parser.add_argument("--topic2", type=str, default="topic2")
    parser.add_argument("--bootstrap_server", type=str, default="127.0.0.1:9092")
    args = parser.parse_args()

    producer1 = threading.Thread(
        target=producer_job,
        args=(args.topic1, stream_1_generator, args.bootstrap_server),
    )
    producer2 = threading.Thread(
        target=producer_job,
        args=(args.topic2, stream_2_generator, args.bootstrap_server),
    )

    producer1.start()
    producer2.start()

    try:
        while True:
            time.sleep(1)
    except KeyboardInterrupt:
        pass
    finally:
        STOP_EVENT.set()
        producer1.join()
        producer2.join()


if __name__ == "__main__":
    main()
