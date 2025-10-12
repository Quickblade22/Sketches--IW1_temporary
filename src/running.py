#!/usr/bin/env python3
"""
Script to run rom_planner C++ executable with binary files location and parameters
"""

import os
import sys
import subprocess
import argparse
import time
import stat
import signal
import threading
import psutil
from datetime import datetime

class RomPlannerRunner:
    def __init__(self, output_file=None):
        self.process = None
        self.shutdown = False
        self.output_file = output_file
        self.child_processes = []
        
    def signal_handler(self, sig, frame):
        """Handle Ctrl+C and other termination signals"""
        print(f"\nReceived signal {sig}, shutting down...")
        self.shutdown = True
        self.terminate_child_processes()
        
        if self.process and self.process.poll() is None:
            print("Terminating rom_planner process...")
            self.process.terminate()
            
            # Wait a bit for graceful termination
            for i in range(5):
                if self.process.poll() is not None:
                    break
                time.sleep(0.5)
            
            # Force kill if still running
            if self.process.poll() is None:
                print("Forcing kill of rom_planner process...")
                self.process.kill()
    
    def setup_signal_handlers(self):
        """Setup signal handlers for graceful shutdown"""
        signal.signal(signal.SIGINT, self.signal_handler)
        signal.signal(signal.SIGTERM, self.signal_handler)
    
    def write_output_to_file(self, line):
        """Write output to file with timestamp"""
        if line and not self.shutdown and self.output_file:
            timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            try:
                with open(self.output_file, 'a', encoding='utf-8') as f:
                    f.write(f"[{timestamp}] {line}\n")
            except Exception as e:
                print(f"Error writing to output file: {e}")

    def terminate_child_processes(self):
        """Terminate all child processes including Arcade Learning Environment"""
        if not self.process:
            return
            
        try:
            # Get the main process
            main_process = psutil.Process(self.process.pid)
            
            # Get all child processes
            children = main_process.children(recursive=True)
            
            # Also look for common ALE processes
            ale_processes = []
            for proc in psutil.process_iter(['pid', 'name', 'cmdline']):
                try:
                    if any(name in proc.info.get('name', '').lower() for name in ['ale', 'arcade', 'emulator']):
                        ale_processes.append(proc)
                    elif proc.info.get('cmdline'):
                        cmdline = ' '.join(proc.info['cmdline']).lower()
                        if any(keyword in cmdline for keyword in ['ale', 'arcade', 'atari', 'emulator']):
                            ale_processes.append(proc)
                except (psutil.NoSuchProcess, psutil.AccessDenied):
                    continue
            
            all_children = children + ale_processes
            
            if all_children:
                print(f"Terminating {len(all_children)} child processes...")
                
                # First try to terminate gracefully
                for child in all_children:
                    try:
                        child.terminate()
                    except (psutil.NoSuchProcess, psutil.AccessDenied):
                        pass
                
                # Wait a bit
                gone, alive = psutil.wait_procs(all_children, timeout=3)
                
                # Force kill any remaining
                for child in alive:
                    try:
                        print(f"Force killing child process {child.pid}...")
                        child.kill()
                    except (psutil.NoSuchProcess, psutil.AccessDenied):
                        pass
                        
        except (psutil.NoSuchProcess, psutil.AccessDenied) as e:
            print(f"Warning: Could not terminate child processes: {e}")
        except Exception as e:
            print(f"Error terminating child processes: {e}")

def fix_executable_permissions(executable_path):
    """Fix executable permissions if needed"""
    if os.path.exists(executable_path):
        # Check if file is executable
        if not os.access(executable_path, os.X_OK):
            print(f"Making {executable_path} executable...")
            try:
                st = os.stat(executable_path)
                os.chmod(executable_path, st.st_mode | stat.S_IEXEC)
                return True
            except Exception as e:
                print(f"Warning: Could not set executable permissions: {e}")
                return False
    return True

def find_rom_planner():
    """Find the rom_planner executable"""
    locations = [
        './rom_planner',
        './src/rom_planner',
        '../src/rom_planner'
    ]
    
    for location in locations:
        if os.path.isfile(location):
            return location
    return None

def run_rom_planner(bin_location, time_budget, simulation_budget, additional_args=None, output_file=None):
    """
    Run the rom_planner executable with specified parameters
    
    Args:
        bin_location (str): Path to binary files directory
        time_budget (float): Time budget parameter
        simulation_budget (int): Simulation budget parameter
        additional_args (list): Additional command line arguments
        output_file (str): Output file for redirection (if None, output goes to console)
    """
    
    runner = RomPlannerRunner(output_file)
    runner.setup_signal_handlers()
    
    # Find the executable
    executable_path = find_rom_planner()
    
    if not executable_path:
        print("Error: rom_planner executable not found!")
        print("Please make sure to compile the project first with 'make'")
        return False
    
    # Fix permissions if needed
    if not fix_executable_permissions(executable_path):
        print(f"Warning: {executable_path} may not be executable")
    
    # Verify the binary files location exists
    if not os.path.exists(bin_location):
        print(f"Error: Binary files location '{bin_location}' does not exist!")
        return False
    
    # Build the command - FIXED THE TYPO: '----simulator-budget' should be '--simulator-budget'
    cmd = [
        executable_path,
        bin_location,
        '--time-budget', str(time_budget),
        '--simulator-budget', str(simulation_budget)  # Fixed: two dashes, not four
    ]
    
    # Add any additional arguments
    if additional_args:
        cmd.extend(additional_args)
    
    print("=" * 60)
    print("Running rom_planner with parameters:")
    print(f"  Binary location: {bin_location}")
    print(f"  Time budget: {time_budget}")
    print(f"  Simulation budget: {simulation_budget}")
    print(f"  Executable: {executable_path}")
    
    # Record start time
    start_time = time.time()
    
    try:
        if output_file:
            # Use shell redirection (>>) to append to output file
            print(f"  Output file: {output_file} (using shell redirection >>)")
            print("=" * 60)
            print("Press Ctrl+C to terminate rom_planner and Arcade Learning Environment")
            print("=" * 60)
            
            # Build the shell command with redirection
            shell_cmd = ' '.join(cmd) + f' >> {output_file} 2>&1'
            
            # Write header to output file
            try:
                with open(output_file, 'a', encoding='utf-8') as f:
                    f.write("=" * 60 + "\n")
                    f.write("rom_planner Execution Log\n")
                    f.write(f"Started at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                    f.write(f"Binary location: {bin_location}\n")
                    f.write(f"Time budget: {time_budget}\n")
                    f.write(f"Simulation budget: {simulation_budget}\n")
                    f.write("=" * 60 + "\n\n")
            except Exception as e:
                print(f"Error writing header to output file: {e}")
            
            print(f"\n--- rom_planner output is being written to {output_file} ---")
            print("The C++ executable output will not appear in the console.")
            print("Check the output file for progress and results.\n")
            
            # Run with shell redirection
            runner.process = subprocess.Popen(
                shell_cmd,
                shell=True,  # Use shell for redirection
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                bufsize=1
            )
            
            # Wait for process to complete or be terminated
            while runner.process.poll() is None and not runner.shutdown:
                time.sleep(0.1)
            
        else:
            # Original behavior: capture output in Python and write to console
            print("  Output: Console")
            print("=" * 60)
            print("Press Ctrl+C to terminate rom_planner and Arcade Learning Environment")
            print("=" * 60)
            
            runner.process = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                bufsize=1,
                universal_newlines=True
            )
            
            # Start threads to read stdout and stderr
            def print_stdout(line):
                if not runner.shutdown:
                    print(line)
            
            def print_stderr(line):
                if not runner.shutdown:
                    print(f"STDERR: {line}", file=sys.stderr)
            
            stdout_thread = threading.Thread(
                target=runner.read_output, 
                args=(runner.process.stdout, print_stdout)
            )
            stderr_thread = threading.Thread(
                target=runner.read_output, 
                args=(runner.process.stderr, print_stderr)
            )
            
            stdout_thread.daemon = True
            stderr_thread.daemon = True
            
            stdout_thread.start()
            stderr_thread.start()
            
            print("\n--- rom_planner Output ---")
            
            # Wait for process to complete or be terminated
            while runner.process.poll() is None and not runner.shutdown:
                time.sleep(0.1)
        
        # If we're shutting down, the signal handler will take care of termination
        if runner.shutdown:
            # Wait a bit more for the process to terminate
            for i in range(10):
                if runner.process.poll() is not None:
                    break
                time.sleep(0.5)
        
        # Get final return code
        return_code = runner.process.poll()
        
        # Calculate execution time
        execution_time = time.time() - start_time
        
        # Write footer to output file if using file output
        if output_file and not runner.shutdown:
            try:
                with open(output_file, 'a', encoding='utf-8') as f:
                    f.write("\n" + "=" * 60 + "\n")
                    f.write("Execution completed!\n")
                    f.write(f"Return code: {return_code}\n")
                    f.write(f"Execution time: {execution_time:.2f} seconds\n")
                    f.write(f"Finished at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                    f.write("=" * 60 + "\n")
            except Exception as e:
                print(f"Error writing footer to output file: {e}")
        
        print("\n" + "=" * 60)
        if runner.shutdown:
            print("Execution terminated by user!")
        else:
            print("Execution completed!")
        print(f"Return code: {return_code}")
        print(f"Execution time: {execution_time:.2f} seconds")
        if output_file:
            print(f"All output has been written to: {output_file}")
        print("=" * 60)
        
        return return_code == 0
        
    except FileNotFoundError:
        print(f"Error: rom_planner executable not found at '{executable_path}'")
        return False
    except Exception as e:
        print(f"Error running rom_planner: {e}")
        return False
    finally:
        # Ensure process is terminated
        if runner.process and runner.process.poll() is None:
            runner.process.terminate()
            runner.process.wait()

    def read_output(self, pipe, output_callback):
        """Read output from a pipe in a separate thread"""
        try:
            with pipe:
                for line in iter(pipe.readline, ''):
                    if self.shutdown:
                        break
                    if line and output_callback:
                        output_callback(line.strip())
        except Exception as e:
            if not self.shutdown:
                print(f"Error reading output: {e}")

def main():
    parser = argparse.ArgumentParser(
        description='Run rom_planner C++ executable with specified parameters',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python run_rom_planner.py /path/to/bin 300.5 10000
  python run_rom_planner.py /data/binaries 600.0 50000 --additional-args="--verbose --output=results.txt"
  python run_rom_planner.py /work/bin 120.0 2000 -a "--threads=4 --mode=fast"
  python run_rom_planner.py /path/to/bin 300.5 10000 --output-file=output.txt
        """
    )
    
    parser.add_argument(
        'bin_location',
        help='Path to binary files directory (required)'
    )
    
    parser.add_argument(
        'time_budget',
        type=float,
        help='Time budget parameter (float, e.g., 300.5 for 300.5 seconds)'
    )
    
    parser.add_argument(
        'simulation_budget',
        type=int,
        help='Simulation budget parameter (integer, e.g., 10000 for 10000 simulations)'
    )
    
    parser.add_argument(
        '-a', '--additional-args',
        help='Additional arguments to pass to rom_planner (as a string)'
    )
    
    parser.add_argument(
        '--compile-first',
        action='store_true',
        help='Compile the project first by running "make"'
    )
    
    parser.add_argument(
        '--output-file',
        help='Output file for shell redirection (uses >> to append)'
    )
    
    args = parser.parse_args()
    
    # Compile first if requested
    if args.compile_first:
        print("Compiling project...")
        try:
            # Use just 'make' without -C
            result = subprocess.run(
                ['make'],
                capture_output=True,
                text=True
            )
            if result.returncode != 0:
                print(f"Compilation failed: {result.stderr}")
                return 1
            print("Compilation successful!")
            
            # Fix permissions after compilation
            rom_planner_path = find_rom_planner()
            if rom_planner_path:
                fix_executable_permissions(rom_planner_path)
                
        except Exception as e:
            print(f"Error during compilation: {e}")
            return 1
    
    # Parse additional arguments
    additional_args = []
    if args.additional_args:
        # Simple split - for more complex cases, consider shlex.split()
        additional_args = args.additional_args.split()
    
    # Run rom_planner
    success = run_rom_planner(
        bin_location=args.bin_location,
        time_budget=args.time_budget,
        simulation_budget=args.simulation_budget,
        additional_args=additional_args,
        output_file=args.output_file
    )
    
    return 0 if success else 1

if __name__ == "__main__":
    sys.exit(main())