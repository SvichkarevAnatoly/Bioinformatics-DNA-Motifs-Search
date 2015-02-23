import argparse


def create_parser():
    parser = argparse.ArgumentParser(description="Matching position weight matrices (PWM) against DNA sequences")
    parser.add_argument("fasta", type=argparse.FileType('r'), help="fasta file with DNA sequences")
    parser.add_argument("pwm", type=argparse.FileType('r'), help="file with position weight matrices (PWM)")
    parser.add_argument("-o", "--output", dest="matching", type=argparse.FileType('w'), metavar='matching',
                        help="output file with matching results")
    return parser


def workflow(args):
    # TODO:
    pass


def main():
    parser = create_parser()
    args = parser.parse_args()
    workflow(args)


if __name__ == "__main__":
    main()