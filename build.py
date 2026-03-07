#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2025 Stephan Willerich
# SPDX-License-Identifier: MIT License

"""
Build script for magnetix_toolbox
Usage: python build.py [options]
Options:
  -j, --jobs <num>    Number of parallel jobs (default: auto-detect)
  -c, --clean         Clean build directory before building
  -v, --verbose       Verbose output
  -h, --help          Show this help message
"""

import os
import sys
import subprocess
import argparse
import shutil
from pathlib import Path
from multiprocessing import cpu_count


class Colors:
    """ANSI color codes for terminal output"""
    BLUE = '\033[0;34m'
    GREEN = '\033[0;32m'
    YELLOW = '\033[1;33m'
    RED = '\033[0;31m'
    NC = '\033[0m'  # No Color

    @staticmethod
    def disable_on_windows():
        """Disable colors on Windows"""
        if sys.platform.startswith('win'):
            for attr in dir(Colors):
                if not attr.startswith('_'):
                    setattr(Colors, attr, '')


class BuildSystem:
    """Orchestrates the multi-component build process"""

    def __init__(self, jobs=None, clean=False, verbose=False):
        self.script_dir = Path(__file__).parent.absolute()
        self.jobs = jobs or cpu_count()
        self.clean = clean
        self.verbose = verbose
        self.install_prefix = self.script_dir / "install"
        self.build_dir = self.script_dir / "build"
        self.cmake_prefix_path = str(self.install_prefix)

    def print_header(self, text):
        """Print a formatted header"""
        print(f"\n{Colors.BLUE}{'='*40}{Colors.NC}")
        print(f"{Colors.BLUE}  {text}{Colors.NC}")
        print(f"{Colors.BLUE}{'='*40}{Colors.NC}\n")

    def print_step(self, step_num, total_steps, text):
        """Print a build step"""
        print(f"{Colors.YELLOW}Step {step_num}/{total_steps}: {text}...{Colors.NC}")

    def print_success(self, text):
        """Print a success message"""
        print(f"{Colors.GREEN}✓ {text}{Colors.NC}")

    def print_error(self, text):
        """Print an error message and exit"""
        print(f"{Colors.RED}✗ Error: {text}{Colors.NC}", file=sys.stderr)
        sys.exit(1)

    def print_warning(self, text):
        """Print a warning message"""
        print(f"{Colors.YELLOW}⚠ {text}{Colors.NC}")

    def check_prerequisites(self):
        """Verify that required tools are available"""
        print(f"{Colors.YELLOW}Checking prerequisites...{Colors.NC}")
        
        required_tools = ['cmake', 'make']
        for tool in required_tools:
            if shutil.which(tool) is None:
                self.print_error(f"{tool} not found in PATH")
        
        self.print_success("Prerequisites OK")
        print()

    def check_dependencies(self):
        """Inform user about dependency requirements"""
        print(f"{Colors.YELLOW}Checking for FEniCSx and dependencies...{Colors.NC}")
        print(f"{Colors.YELLOW}(These should be in your environment or CMAKE_PREFIX_PATH){Colors.NC}")
        
        if shutil.which('module') is not None:
            print(f"{Colors.YELLOW}Tip: If not available, try: module load dolfinx{Colors.NC}")
        
        print()

    def run_command(self, cmd, cwd=None, description=None):
        """Execute a shell command"""
        if self.verbose:
            print(f"{Colors.BLUE}Running: {' '.join(cmd)}{Colors.NC}")
        
        try:
            output = subprocess.PIPE if not self.verbose else None
            result = subprocess.run(
                cmd,
                cwd=cwd,
                stdout=output,
                stderr=subprocess.STDOUT,
                check=True,
                text=True
            )
            if self.verbose and result.stdout:
                print(result.stdout)
            return True
        except subprocess.CalledProcessError as e:
            if description:
                self.print_error(f"{description}\n{e.output if e.output else ''}")
            else:
                self.print_error(f"Command failed: {' '.join(cmd)}\n{e.output if e.output else ''}")
            return False

    def clean_build_dirs(self):
        """Remove build and install directories"""
        print(f"{Colors.YELLOW}Cleaning build directories...{Colors.NC}")
        
        for d in [self.build_dir, self.install_prefix]:
            if d.exists():
                shutil.rmtree(d)
        
        self.print_success("Build directories cleaned")

    def build_component(self, step_num, total_steps, component_name, source_dir):
        """Build a single component"""
        self.print_step(step_num, total_steps, f"Building {component_name}")
        
        build_subdir = self.build_dir / f"{component_name}_build"
        build_subdir.mkdir(parents=True, exist_ok=True)
        
        # Configure
        cmake_cmd = [
            'cmake',
            f'-DCMAKE_INSTALL_PREFIX={self.install_prefix}',
            f'-DCMAKE_PREFIX_PATH={self.cmake_prefix_path}',
            str(source_dir)
        ]
        if self.verbose:
            cmake_cmd.insert(1, '--verbose')
        
        self.run_command(cmake_cmd, cwd=build_subdir, description=f"CMake configuration for {component_name}")
        
        # Build
        build_cmd = ['cmake', '--build', '.', '-j', str(self.jobs)]
        if self.verbose:
            build_cmd.append('--verbose')
        
        self.run_command(build_cmd, cwd=build_subdir, description=f"Build failed for {component_name}")
        
        # Install (if not examples)
        if component_name != 'examples':
            install_cmd = ['cmake', '--build', '.', '--target', 'install']
            if self.verbose:
                install_cmd.append('--verbose')
            
            self.run_command(install_cmd, cwd=build_subdir, description=f"Install failed for {component_name}")
            
            # Update prefix path for next components
            if component_name == 'hysteresis_model':
                self.cmake_prefix_path = f"{self.install_prefix}:{self.cmake_prefix_path}"
        
        self.print_success(f"{component_name} built successfully")
        print()
        
        # Return build directory for examples
        if component_name == 'examples':
            return build_subdir

    def organize_artifacts(self, examples_build_dir):
        """Place built artifacts in convenient locations"""
        print(f"{Colors.YELLOW}Organizing build artifacts...{Colors.NC}")
        
        exec_src = examples_build_dir / "magnetostatic_2D_exec"
        exec_dest = self.script_dir / "examples" / "magnetostatic_2D" / "magnetostatic_2D_exec"
        
        if exec_src.exists():
            shutil.copy2(exec_src, exec_dest)
            self.print_success(f"Executable placed at {exec_dest.relative_to(self.script_dir)}")
        
        print()

    def print_environment_setup(self):
        """Print environment setup instructions"""
        hysteresis_path = self.install_prefix / "dpc_hysteresis-0.1"
        fenicsx_path = self.install_prefix / "fenicsx_magnetics_toolbox-0.10"
        
        print(f"{Colors.YELLOW}To setup your environment, add to ~/.bashrc or ~/.zshrc:{Colors.NC}\n")
        print(f"{Colors.GREEN}export CMAKE_PREFIX_PATH=\"{hysteresis_path}:{fenicsx_path}:$CMAKE_PREFIX_PATH\"{Colors.NC}")
        print(f"{Colors.GREEN}export PYTHONPATH=\"{fenicsx_path}/python:$PYTHONPATH\"{Colors.NC}")
        
        print(f"\n{Colors.YELLOW}Or run to update your current shell:{Colors.NC}\n")
        print(f"{Colors.GREEN}export CMAKE_PREFIX_PATH=\"{hysteresis_path}:{fenicsx_path}:$CMAKE_PREFIX_PATH\"{Colors.NC}")
        print(f"{Colors.GREEN}export PYTHONPATH=\"{fenicsx_path}/python:$PYTHONPATH\"{Colors.NC}")
        
        print(f"\n{Colors.YELLOW}Run TEAM Problem 32 example:{Colors.NC}\n")
        print(f"{Colors.GREEN}cd examples/magnetostatic_2D/TEAM_Problem_32{Colors.NC}")
        print(f"{Colors.GREEN}mpirun -np 4 ../magnetostatic_2D_exec --scen TeamProblem32_case3_default.xml{Colors.NC}")
        print()

    def build(self):
        """Execute the complete build process"""
        self.print_header("magnetix_toolbox Build Script")
        
        self.check_prerequisites()
        self.check_dependencies()
        
        if self.clean:
            self.clean_build_dirs()
        
        # Create build directory
        if not self.build_dir.exists():
            print(f"{Colors.YELLOW}Creating build directory...{Colors.NC}")
            self.build_dir.mkdir(parents=True)
            self.print_success("Build directory created")
            print()
        
        # Build components
        examples_dir = self.build_component(
            1, 3, 'hysteresis_model',
            self.script_dir / 'hysteresis_model'
        )
        
        self.build_component(
            2, 3, 'fenicsx_tools',
            self.script_dir / 'fenicsx_tools' / 'library'
        )
        
        examples_build_dir = self.build_component(
            3, 3, 'examples',
            self.script_dir / 'examples' / 'magnetostatic_2D'
        )
        
        self.organize_artifacts(examples_build_dir)
        
        self.print_header("Build successful!")
        self.print_environment_setup()


def main():
    """Main entry point"""
    # Disable colors on Windows
    Colors.disable_on_windows()
    
    parser = argparse.ArgumentParser(
        description='Build script for magnetix_toolbox',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    
    parser.add_argument(
        '-j', '--jobs',
        type=int,
        default=None,
        help='Number of parallel jobs (default: auto-detect)'
    )
    parser.add_argument(
        '-c', '--clean',
        action='store_true',
        help='Clean build directory before building'
    )
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Verbose output'
    )
    
    args = parser.parse_args()
    
    builder = BuildSystem(
        jobs=args.jobs,
        clean=args.clean,
        verbose=args.verbose
    )
    
    try:
        builder.build()
    except KeyboardInterrupt:
        print(f"\n{Colors.RED}Build interrupted by user{Colors.NC}", file=sys.stderr)
        sys.exit(130)
    except Exception as e:
        print(f"{Colors.RED}Unexpected error: {e}{Colors.NC}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
