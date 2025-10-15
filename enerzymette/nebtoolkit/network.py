from typing import Optional
import socket


def find_available_port(start_port: int = 5000, max_port: int = 65535, host: str = "0.0.0.0") -> Optional[int]:
    """
    Find an available port starting from a given port number.
    
    Args:
        start_port: The port number to start searching from (default: 5000)
        max_port: The maximum port number to check (default: 65535)
        host: The host address to bind to (default: "0.0.0.0")
    
    Returns:
        The first available port number, or None if no port is available
    
    Example:
        >>> port = find_available_port()
        >>> print(f"Available port: {port}")
        Available port: 5000
        
        >>> port = find_available_port(start_port=8000)
        >>> print(f"Available port: {port}")
        Available port: 8000
    """
    for port in range(start_port, max_port + 1):
        try:
            # Try to create a socket and bind to the port
            with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
                # Set socket option to reuse address
                sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
                # Try to bind to the port
                sock.bind((host, port))
                # If successful, the port is available
                return port
        except OSError:
            # Port is in use or cannot be bound, try the next one
            continue
    
    # No available port found
    return None