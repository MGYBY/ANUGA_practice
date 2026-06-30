"""Run the two-floor Kim et al. ATM demonstration in sequence."""
from run_ground import run_ground
from run_basement import run_basement


def main() -> None:
    run_ground()
    run_basement()


if __name__ == "__main__":
    main()
