if [[ ":$PYTHONPATH:" != *":/hppfs/work/pn69ni/di73jef3/Softwares/Kratos/bin/Release:"* ]]; then
    export PYTHONPATH=/hppfs/work/pn69ni/di73jef3/Softwares/Kratos/bin/Release:$PYTHONPATH
    export LD_LIBRARY_PATH=/hppfs/work/pn69ni/di73jef3/Softwares/Kratos/bin/Release/libs:$LD_LIBRARY_PATH
    echo "adding /hppfs/work/pn69ni/di73jef3/Softwares/Kratos/bin/Release to PYTHONPATH"
    echo "adding /hppfs/work/pn69ni/di73jef3/Softwares/Kratos/bin/Release/libs to LD_LIBRARY_PATH"
fi
echo "adding paths finished"