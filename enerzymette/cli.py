import argparse


def get_parser():
    parser = argparse.ArgumentParser(
        description="Enerzymette cli",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    subparsers = parser.add_subparsers(title="Valid subcommands", dest="command")

    parser_idpp = subparsers.add_parser(
        "idpp",
        help="IDPP interpolation between two molecules",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser_idpp.add_argument('-r', '--reactant', type=str, 
        help='input reactant xyz file path'
    )
    parser_idpp.add_argument('-p', '--product', type=str, 
        help='input product xyz file path'
    )
    parser_idpp.add_argument('-o', '--output', type=str, 
        help='output interpolated xyz file path'
    )
    parser_idpp.add_argument('-n', '--n_images', type=int,
        help='number of images to generate'
    )
    parser_idpp.add_argument('-c', '--constraints', type=str,
        help='terachem input file path'
    )
    
    parser_terachem_timing = subparsers.add_parser(
        "terachem_timing",
        help="Calculate the wall time of a terachem job",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser_terachem_timing.add_argument('-f', '--filename', type=str,
        help='terachem output file path'
    )

    parser_orca_terachem_request = subparsers.add_parser(
        "orca_terachem_request",
        help="Request a terachem job from orca",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser_orca_terachem_request.add_argument('-i', '--input', type=str,
        help='orca input file path'
    )
    parser_orca_terachem_request.add_argument('-t', '--template', type=str,
        help='terachem input template file path'
    )

    parser_launch_enerzyme_neb = subparsers.add_parser(
        "enerzyme_neb",
        help="Launch a enerzyme neb job",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser_launch_enerzyme_neb.add_argument('-r', '--reactant', type=str,
        help='initial reactant path'
    )
    parser_launch_enerzyme_neb.add_argument('-p', '--product', type=str,
        help='initial product path'
    )
    parser_launch_enerzyme_neb.add_argument('-o', '--output', type=str,
        help='output path', default="."
    )
    parser_launch_enerzyme_neb.add_argument('-m', '--model', type=str,
        help='model path', default=".."
    )
    parser_launch_enerzyme_neb.add_argument('-q', '--reference', type=str,
        help='quantum chemistry parameters reference path'
    )
    parser_launch_enerzyme_neb.add_argument('-c', '--server_config', type=str,
        help='server config path'
    )
    parser_launch_enerzyme_neb.add_argument('-n', '--n_images', type=int,
        help='number of images', default=25
    )
    parser_launch_enerzyme_neb.add_argument('-b', '--port', type=int,
        help='port', default=5000
    )
    parser_launch_enerzyme_neb.add_argument('-i', '--interrupt_strategy', type=str,
        help='interrupt strategy', default="stdout"
    )

    args = parser.parse_args()
    return args


def main():
    args = get_parser()
    if args.command == 'idpp':
        from .idpp import idpp
        return idpp(
            reactant_path=args.reactant,
            product_path=args.product,
            output_path=args.output,
            n_images=args.n_images,
            constraints_file=args.constraints,
        )
    elif args.command == 'terachem_timing':
        from .terachem.timing import terachem_timing
        return terachem_timing(
            filename=args.filename
        )
    elif args.command == 'orca_terachem_request':
        from .terachem.orca_io import run
        return run(
            orca_extinp_file=args.input,
            terachem_input_template=args.template,
        )
    elif args.command == "enerzyme_neb":
        from .nebtoolkit.launcher import EnerzymeNEBLauncher
        launcher = EnerzymeNEBLauncher(
            reactant_path=args.reactant,
            product_path=args.product,
            output_path=args.output,
            model_path=args.model,
            reference_path=args.reference,
            server_config_path=args.server_config,
            n_images=args.n_images,
            port=args.port,
            interrupt_strategy=args.interrupt_strategy,
        )
        launcher.launch()
    else:
        raise NotImplementedError(f"Command {args.command} is not supported now.")


if __name__ == '__main__':
    main()