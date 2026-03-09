import shutil
import os
import socket

def check_disk_usage(path):
    """Checks the disk usage of a given path."""
    print(f"🔍 Checking path: {path}")
    if os.path.exists(path):
        usage = shutil.disk_usage(path)
        print(f"📊 Disk Status for {path} on {socket.gethostname()}:")
        print(f"   Total: {usage.total / (1024**3):.2f} GB")
        print(f"   Used:  {usage.used / (1024**3):.2f} GB")
        print(f"   Free:  {usage.free / (1024**3):.2f} GB")
        
        # Check if the directory is actually mounted (Linux specific check)
        is_mount = os.path.ismount(path)
        print(f"🔗 Is active mount point: {'✅ YES' if is_mount else '❌ NO'}")
        
        # Test write permissions
        test_file = os.path.join(path, f".test_write_{socket.gethostname()}")
        try:
            with open(test_file, 'w') as f:
                f.write('test')
            os.remove(test_file)
            print("✍️  Write Permissions: ✅ OK")
        except Exception as e:
            print(f"✍️  Write Permissions: ❌ FAILED ({e})")
    else:
        print(f"⚠️  Path {path} NOT FOUND on this node!")

if __name__ == "__main__":
    # Check both the local/master mount and the potential worker mount
    paths_to_check = ["/media/nfs01", "/mnt/astro_data"]
    
    print("--- 📡 Storage Connectivity Check ---")
    for p in paths_to_check:
        check_disk_usage(p)
        print("-" * 38)
