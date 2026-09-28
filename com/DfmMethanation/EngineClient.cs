using System;
using System.Diagnostics;
using System.IO;
using System.Text;

namespace PhD.DfmMethanation
{
    /// <summary>
    /// Runs the step-8 Docker engine: cape_request.json → cape_result.json.
    /// FEniCSx stays out of the PME process.
    /// </summary>
    internal static class EngineClient
    {
        public static string Run(string workDir, string image, string extraArgs)
        {
            var req = Path.Combine(workDir, "cape_request.json");
            var res = Path.Combine(workDir, "cape_result.json");
            if (!File.Exists(req))
                throw new InvalidOperationException("cape_request.json missing");
            if (File.Exists(res))
                File.Delete(res);

            var args = new StringBuilder();
            args.Append("run --rm ");
            args.Append("-v \"").Append(workDir).Append(":/data\" ");
            args.Append(image).Append(' ');
            args.Append("/data/cape_request.json -o /data/cape_result.json");
            if (!string.IsNullOrWhiteSpace(extraArgs))
                args.Append(' ').Append(extraArgs);

            var psi = new ProcessStartInfo
            {
                FileName = "docker",
                Arguments = args.ToString(),
                UseShellExecute = false,
                RedirectStandardOutput = true,
                RedirectStandardError = true,
                CreateNoWindow = true
            };
            using (var p = Process.Start(psi))
            {
                if (p == null)
                    throw new InvalidOperationException("failed to start docker");
                var stdout = p.StandardOutput.ReadToEnd();
                var stderr = p.StandardError.ReadToEnd();
                p.WaitForExit();
                if (p.ExitCode != 0)
                    throw new InvalidOperationException(
                        "Docker CAPE engine failed: " + stderr + stdout);
            }
            if (!File.Exists(res))
                throw new InvalidOperationException("engine did not write cape_result.json");
            return File.ReadAllText(res, Encoding.UTF8);
        }
    }
}
