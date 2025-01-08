import os
import logging
import webbrowser
from http.server import HTTPServer, SimpleHTTPRequestHandler
from threading import Thread

def create_igv_html(output_dir, out_prefix, ref_genome, target_region):
    template_path = os.path.join(os.path.dirname(__file__), "resources/igv_viewer_template.html")
    
    with open(template_path, "r") as template_file:
        html_content = template_file.read()
    
    logging.info("Creating IGV viewer HTML file")

    # Debugging: Check values before creating HTML
    logging.debug(f"Reference Genome: {ref_genome}")
    logging.debug(f"Target Region: {target_region}")


    # Replace placeholders with actual values
    html_content = html_content.replace("{{ref_genome}}", ref_genome)
    html_content = html_content.replace("{{target_region}}", target_region)
    html_content = html_content.replace("{{bam_file}}", f"{out_prefix}.bam")
    
    output_path = os.path.join(output_dir, "igv_viewer.html")
    with open(output_path, "w") as output_file:
        output_file.write(html_content)

def start_http_server(directory, port=8000):
    os.chdir(directory)
    handler = SimpleHTTPRequestHandler
    httpd = HTTPServer(('localhost', port), handler)
    logging.info(f"Serving HTTP on localhost port {port}")
    httpd.serve_forever()

def open_igv_viewer(output_dir):
    # Start HTTP server in a separate thread
    server_thread = Thread(target=start_http_server, args=(output_dir,))
    server_thread.daemon = True  # Set as a daemon thread
    server_thread.start()
    logging.info(f"Starting IGV viewer on http://localhost:8000")

    # Open the browser
    webbrowser.open(f"http://localhost:8000/igv_viewer.html")

    # Keep the main thread alive while the server is running
    try:
        while server_thread.is_alive():
            server_thread.join(1)  # Wait for the server thread to finish
    except KeyboardInterrupt:
        logging.info("Shutting down server...")