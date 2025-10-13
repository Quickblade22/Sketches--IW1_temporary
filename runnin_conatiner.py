#!/usr/bin/env python3
import subprocess
import sys
import os
import signal
import time
import os.path

class ContainerRunner:
    def __init__(self):
        self.processes = []
        
    def run_command(self, cmd, shell=False):
        """Run a command and return success status"""
        try:
            print(f"Running: {' '.join(cmd) if isinstance(cmd, list) else cmd}")
            result = subprocess.run(cmd, shell=shell, check=True)
            return result.returncode == 0
        except subprocess.CalledProcessError as e:
            print(f"Error running command: {e}")
            return False
        except KeyboardInterrupt:
            print("Command interrupted by user")
            return False
    
    def build_container(self):
        """Build the Apptainer container"""
        print("Building Apptainer container...")
        return self.run_command(["apptainer", "build", "mycontainer.sif", "apptainers.def"])
    
    def run_in_container_exec(self):
        """Run using apptainer exec for better process control"""
        try:
            print("Running Python script in container using exec...")
            
            # Change to the src directory and run the Python script
            cmd = [
                "apptainer", "exec", 
                "--bind", "/root/Sketches--IW1_temporary",
                "mycontainer.sif",
                "sh", "-c",
                "cd /root/Sketches--IW1_temporary/src && python3 running.py /root/Adventure.bin 10000 5000000 --output-file=output.txt"
            ]
            
            proc = subprocess.Popen(cmd)
            self.processes.append(proc)
            
            # Wait for process to complete
            returncode = proc.wait()
            self.processes.remove(proc)
            
            if returncode == 0:
                print("Container execution completed successfully!")
                return True
            else:
                print(f"Container execution failed with return code: {returncode}")
                return False
                
        except Exception as e:
            print(f"Error during container execution: {e}")
            return False
    
    def cleanup(self):
        """Clean up all running processes"""
        for proc in self.processes:
            try:
                proc.terminate()
                proc.wait(timeout=5)
            except:
                try:
                    proc.kill()
                except:
                    pass
    
    def signal_handler(self, sig, frame):
        """Handle termination signals"""
        print(f"\nReceived signal {sig}, cleaning up...")
        self.cleanup()
        sys.exit(0)

def main():
    """Main function"""
    runner = ContainerRunner()
    
    # Set up signal handling
    signal.signal(signal.SIGINT, runner.signal_handler)
    signal.signal(signal.SIGTERM, runner.signal_handler)
    
    print("Starting Apptainer build and execution script...")
    if(not os.path.isfile("mycontainer.sif")):
        # Step 1: Build the container
        if not runner.build_container():
            print("Container build failed. Exiting.")
            sys.exit(1)
        
        print("Container built successfully!")
    
    # Step 2: Run in the container using exec (recommended for better control)
    if not runner.run_in_container_exec():
        print("Container execution failed. Exiting.")
        sys.exit(1)
    
    print("Script completed successfully!")

if __name__ == "__main__":
    main()