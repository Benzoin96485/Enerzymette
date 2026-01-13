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
    parser_launch_enerzyme_neb.add_argument('--optimization_method', type=str,
        help='optimization method', default="LBFGS"
    )
    parser_launch_enerzyme_neb.add_argument('-i', '--interrupt_strategy', type=str,
        help='interrupt strategy', default="stdout"
    )
    parser_launch_enerzyme_neb.add_argument('--max_restart_attempts', type=int,
        help='max restart attempts', default=10
    )
    parser_launch_enerzyme_neb.add_argument('--min_spring_constant', type=float,
        help='min spring constant', default=0.01
    )
    parser_launch_enerzyme_neb.add_argument('--max_spring_constant', type=float,
        help='max spring constant', default=0.1
    )

    parser_launch_enerzyme_scan = subparsers.add_parser(
        "enerzyme_scan",
        help="Launch a enerzyme scan job",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser_launch_enerzyme_scan.add_argument('-r', '--reactant', type=str,
        help='initial reactant path'
    )
    parser_launch_enerzyme_scan.add_argument('-o', '--output', type=str,
        help='output path', default="."
    )
    parser_launch_enerzyme_scan.add_argument('-m', '--model', type=str,
        help='model path', default=".."
    )
    parser_launch_enerzyme_scan.add_argument('-q', '--reference', type=str,
        help='quantum chemistry parameters reference path'
    )
    parser_launch_enerzyme_scan.add_argument('-n', '--n_steps', type=int,
        help='number of steps', default=25
    )
    parser_update_terachem_scan = subparsers.add_parser(
        "update_terachem_scan",
        help="Update a terachem scan input file",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser_update_terachem_scan.add_argument('-i', '--terachem_input', type=str,
        help='terachem input file path'
    )
    parser_update_terachem_scan.add_argument('-s', '--updated_structure', type=str,
        help='updated structure xyz file path'
    )
    parser_update_terachem_scan.add_argument('-o', '--output_path', type=str,
        help='output path', default="."
    )

    parser_enerzyme_active_learning = subparsers.add_parser(
        "enerzyme_active_learning",
        help="Launch a enerzyme active learning job",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser_enerzyme_active_learning.add_argument('-p', '--pretrain', type=str,
        help='pretrain path'
    )
    parser_enerzyme_active_learning.add_argument('-o', '--output', type=str,
        help='output path', default="."
    )
    parser_enerzyme_active_learning.add_argument('-t', '--tmp', type=str,
        help='tmp path', default="."
    )
    parser_enerzyme_active_learning.add_argument('-cp', '--calculator_patch', type=str,
        help='calculator patch'
    )
    parser_enerzyme_active_learning.add_argument('-pp', '--plumed_patch', type=str,
        help='plumed patch'
    )   
    parser_enerzyme_active_learning.add_argument('-sc', '--simulation_config', type=str,
        help='simulation config path'
    )
    parser_enerzyme_active_learning.add_argument('-ec', '--extraction_config', type=str,
        help='extraction config path'
    )
    parser_enerzyme_active_learning.add_argument('-ac', '--annotation_config', type=str,
        help='annotation config path'
    )
    parser_enerzyme_active_learning.add_argument('-tc', '--training_config', type=str,
        help='training config path'
    )
    parser_enerzyme_active_learning.add_argument('-n', '--n_iterations', type=int,
        help='number of iterations', default=10
    )
    parser_enerzyme_active_learning.add_argument('-r', '--training_ratio', type=float,
        help='training ratio', default=0.8
    )
    parser_enerzyme_active_learning.add_argument('-np', '--n_presimulation_steps_per_iteration', type=int,
        help='number of presimulation steps per iteration', default=0
    )
    parser_enerzyme_active_learning.add_argument('-b', '--cluster_inference_batch_size', type=int,
        help='cluster inference batch size', default=4
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
            optimization_method=args.optimization_method,
            max_restart_attempts=args.max_restart_attempts,
            min_spring_constant=args.min_spring_constant,
            max_spring_constant=args.max_spring_constant
        )
        launcher.launch()
    elif args.command == "enerzyme_scan":
        from .scantoolkit.launcher import EnerzymeScanLauncher
        launcher = EnerzymeScanLauncher(
            reactant_path=args.reactant,
            output_path=args.output,
            model_path=args.model,
            reference_path=args.reference,
            n_steps=args.n_steps,
        )
        launcher.launch()
    elif args.command == "update_terachem_scan":
        from .scantoolkit.io import update_terachem_scan_input
        update_terachem_scan_input(
            terachem_input_file=args.terachem_input,
            updated_structure_xyz=args.updated_structure,
            output_path=args.output_path,
        )
    elif args.command == "enerzyme_active_learning":
        from .altoolkit.launcher import active_learning_launcher
        launcher = active_learning_launcher(
            pretrain_path=args.pretrain,
            output_path=args.output,
            tmp_path=args.tmp,
            calculator_patch_key=args.calculator_patch,
            plumed_patch_key=args.plumed_patch,
            simulation_config_path=args.simulation_config,
            extraction_config_path=args.extraction_config,
            annotation_config_path=args.annotation_config,
            training_config_path=args.training_config,
            n_iterations=args.n_iterations,
            training_ratio=args.training_ratio,
            cluster_inference_batch_size=args.cluster_inference_batch_size,
            n_presimulation_steps_per_iteration=args.n_presimulation_steps_per_iteration,
        )
        launcher.launch()
    else:
        raise NotImplementedError(f"Command {args.command} is not supported now.")


if __name__ == '__main__':
    main()