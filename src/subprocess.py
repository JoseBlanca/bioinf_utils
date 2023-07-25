import tempfile
import subprocess


def run_in_sh_script(cmd, cwd=None, capture_output=False):
    with tempfile.NamedTemporaryFile(suffix=".sh", mode="tw") as shell_fhand:
        cmd_str = " ".join(map(str, cmd))
        shell_fhand.write(cmd_str)
        shell_fhand.write("\nexit $?")
        shell_fhand.flush()
        try:
            subprocess.run(
                ["sh", shell_fhand.name],
                check=True,
                cwd=cwd,
                capture_output=capture_output,
            )
        except subprocess.CalledProcessError:
            print(" ".join(cmd))
            raise
