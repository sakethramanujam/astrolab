import dask.array as da
from distributed import Client
import time

def process_astro_task(x):
    # Imagine a complex image reduction task
    time.sleep(1) # Simulate heavy processing
    return x ** 2

if __name__ == "__main__":
    # --- 1. Connect to Local/Network Cluster ---
    # To use your network nodes, run 'dask-scheduler' and 'dask-worker' on them.
    # Replace 'localhost' with your master IP if running across network.
    try:
        client = Client('127.0.0.1:8786', timeout="2s")
        print("🔗 Connected to existing Dask cluster on network!")
    except:
        client = Client() # Falls back to local cores
        print("🏠 No network cluster found. Falling back to LOCAL cores.")
    
    print(f"📊 Cluster Info: {client.dashboard_link}")

    # --- 2. Parallelizing a Large Processing Job ---
    # Imagine 1000 FITS images to process
    data = da.random.random((1000, 1000, 100), chunks=(100, 100, 10))
    
    print("🚀 Parallelizing astronomical pixel math across the network...")
    result = data.mean(axis=2).compute()
    
    print(f"✅ Finished! Resulting image shape: {result.shape}")
    client.close()
