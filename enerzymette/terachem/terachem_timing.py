#!/usr/bin/env python3
"""
Comprehensive script to calculate total wall time from TeraChem output files.
Handles multiple iteration tables in a single file.
Marks unfinished tables for on-the-fly analysis.
"""

def is_table_unfinished(lines, start_line, next_table_start):
    """
    Heuristic to determine if a table is unfinished:
    - If it's the last table and the file does not end with a blank line, summary, or 'FINAL ENERGY', it's likely unfinished.
    - If the table is interrupted (not followed by a blank line or summary), it's likely unfinished.
    """
    end_line = next_table_start if next_table_start is not None else len(lines)
    # Look for 'FINAL ENERGY' or 'Converged' or a blank line after the table
    for i in range(end_line-1, start_line, -1):
        line = lines[i].strip()
        if not line:
            return False  # blank line after table, likely finished
        if 'FINAL ENERGY' in line or 'Converged' in line or 'converged' in line:
            return False
        if 'Summary' in line or 'summary' in line:
            return False
    return True  # No sign of completion found

def calculate_wall_time_complete(filename, save_to_file=None):
    """Calculate total wall time from TeraChem output file, handling multiple tables and unfinished tables."""
    try:
        with open(filename, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
        return None
    
    # Find all iteration table headers
    table_starts = []
    for i, line in enumerate(lines):
        if 'Iter' in line and 'DIIS Error' in line and 'Time(s)' in line:
            table_starts.append(i)
    
    # Prepare output
    output_lines = []
    output_lines.append(f"TeraChem Wall Time Analysis: {filename}")
    output_lines.append(f"Found {len(table_starts)} iteration tables")
    output_lines.append("=" * 90)
    
    all_times = []
    total_time = 0.0
    total_iterations = 0
    
    # Print header for table
    output_lines.append(f"{'Table':<6} {'Iterations':<12} {'Wall Time (s)':<15} {'Wall Time (min)':<15} {'Avg/Iter (s)':<12} {'Status':<15}")
    output_lines.append("-" * 90)
    
    for table_idx, start_line in enumerate(table_starts):
        table_times = []
        table_iterations = 0
        next_table_start = table_starts[table_idx+1] if table_idx+1 < len(table_starts) else None
        
        # Process lines after the header until we hit a line that doesn't match iteration pattern
        for i in range(start_line + 1, len(lines)):
            line = lines[i].strip()
            # Stop if we hit another table header or end of file
            if ('Iter' in line and 'DIIS Error' in line) or not line:
                break
            # Look for iteration lines (start with | and have enough columns)
            if line.startswith('|') and len(line.split()) >= 9:
                parts = line.split()
                try:
                    iteration = int(parts[1])
                    time = float(parts[-1])
                    table_times.append((iteration, time))
                    table_iterations = max(table_iterations, iteration)
                except (ValueError, IndexError):
                    continue
        # Heuristic: unfinished if at end of file and no summary, or if interrupted
        unfinished = is_table_unfinished(lines, start_line, next_table_start) if table_times else False
        status = '(unfinished)' if unfinished else ''
        if table_times:
            table_total = sum(time for _, time in table_times)
            total_time += table_total
            total_iterations += table_iterations
            all_times.extend(table_times)
            avg_time = table_total / len(table_times)
            output_lines.append(f"{table_idx + 1:<6} {table_iterations:<12} {table_total:<15.2f} {table_total/60:<15.2f} {avg_time:<12.2f} {status:<15}")
        else:
            output_lines.append(f"{table_idx + 1:<6} {'0':<12} {'0.00':<15} {'0.00':<15} {'0.00':<12} {'':<15}")
    # Print summary
    output_lines.append("-" * 90)
    if all_times:
        overall_avg = total_time / len(all_times)
        times_only = [time for _, time in all_times]
        min_time = min(times_only)
        max_time = max(times_only)
        output_lines.append(f"{'TOTAL':<6} {total_iterations:<12} {total_time:<15.2f} {total_time/60:<15.2f} {overall_avg:<12.2f} {'':<15}")
        output_lines.append("=" * 90)
        output_lines.append(f"Summary Statistics:")
        output_lines.append(f"  • Total iterations across all tables: {total_iterations}")
        output_lines.append(f"  • Total wall time: {total_time:.2f} seconds ({total_time/60:.2f} minutes)")
        output_lines.append(f"  • Average time per iteration: {overall_avg:.2f} seconds")
        output_lines.append(f"  • Time range per iteration: {min_time:.2f} - {max_time:.2f} seconds")
        output_lines.append(f"  • Number of tables processed: {len(table_starts)}")
        # Print to console
        for line in output_lines:
            print(line)
        # Save to file if requested
        if save_to_file:
            try:
                with open(save_to_file, 'w') as f:
                    for line in output_lines:
                        f.write(line + '\n')
                print(f"\nResults saved to: {save_to_file}")
            except Exception as e:
                print(f"Warning: Could not save to file '{save_to_file}': {e}")
        return total_time
    else:
        output_lines.append(f"No iteration data found in {filename}")
        for line in output_lines:
            print(line)
        return None

def terachem_timing(filename):
    calculate_wall_time_complete(filename)
