using System;
using System.Globalization;
using System.IO;
using System.Runtime.InteropServices;

namespace PhD.DfmMethanation
{
    /// <summary>
    /// Append-only trace of COFE calls. Each line is flushed so the last
    /// entry survives an access violation in COFE.exe, which C# cannot catch.
    /// </summary>
    internal static class CapeTrace
    {
        static readonly object Gate = new object();

        public static readonly string LogPath =
            Path.Combine(Path.GetTempPath(), "dfm-cape-open.log");

        static CapeTrace()
        {
            try
            {
                AppDomain.CurrentDomain.UnhandledException += (sender, args) =>
                    Write("UNHANDLED " + args.ExceptionObject);
            }
            catch
            {
                // tracing must not block COM activation
            }
            Write("log " + LogPath);
        }

        public static void Write(string message)
        {
            var line = DateTime.Now.ToString("HH:mm:ss.fff ", CultureInfo.InvariantCulture)
                + message + Environment.NewLine;
            lock (Gate)
            {
                try
                {
                    File.AppendAllText(LogPath, line);
                }
                catch
                {
                    // never throw from the logger
                }
            }
        }

        public static void Guard(string name, Action action)
        {
            Guard(name, () =>
            {
                action();
                return true;
            });
        }

        public static T Guard<T>(string name, Func<T> func)
        {
            Write(">> " + name);
            try
            {
                T value = func();
                Write("<< " + name);
                return value;
            }
            catch (COMException ex)
            {
                Write("!! " + name + " COM 0x" + ex.ErrorCode.ToString("X8") + " " + ex.Message);
                throw;
            }
            catch (Exception ex)
            {
                Write("!! " + name + Environment.NewLine + ex);
                throw new COMException(name + ": " + ex.Message, ex);
            }
        }
    }
}
