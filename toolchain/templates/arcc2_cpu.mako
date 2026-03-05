#!/usr/bin/env bash

<%namespace name="helpers" file="helpers.mako"/>
<%
CORES_MULTIPLIER = 1
%>

% if engine == 'batch':
#SBATCH --nodes=${nodes}
#SBATCH --cpus-per-task=1
#SBATCH --job-name="${name}"
#SBATCH --time=${walltime}
#SBATCH --output="${name}"_%j.out
#SBATCH --error="${name}"_%j.err
#SBATCH --export=ALL
#SBATCH --ntasks-per-node=${tasks_per_node*CORES_MULTIPLIER}
% if partition:
#SBATCH --partition=${partition}
% endif
% if account:
#SBATCH --account="${account}"
% endif
% if quality_of_service:
#SBATCH --qos="${quality_of_service}"
% endif
% if email:
#SBATCH --mail-user=${email}
#SBATCH --mail-type="BEGIN,END,FAIL"
% endif
% endif

${helpers.template_prologue()}

ok ":) Loading modules:\n"
cd "${MFC_ROOT_DIR}"
. ./mfc.sh load -c x -m c
cd - > /dev/null
echo

export LD_LIBRARY_PATH=/opt/ohpc/pub/libs/usr/local/lib64:$LD_LIBRARY_PATH

% for target in targets:
    ${helpers.run_prologue(target)}

    % if not mpi:
        (set -x; ${profiler} "${target.get_install_binpath(case)}")
    % else:
        (set -x; ${profiler} \
            mpirun -np ${nodes*tasks_per_node*CORES_MULTIPLIER} \
                   --bind-to none \
                   "${target.get_install_binpath(case)}")
    % endif

    ${helpers.run_epilogue(target)}

    echo
% endfor

${helpers.template_epilogue()}
