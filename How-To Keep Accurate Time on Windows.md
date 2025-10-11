# How-To: Accurate Time Keeping on Windows 11 for Astronomy

## Introduction
For astronomy applications precise timestamps are essential for tracking fast-moving objects or transient light sources. Windows 11's built-in time synchronization can be optimized for millisecond-level accuracy without additional hardware, by configuring the Windows Time service (W32Time). This guide covers defaults and recommended changes.

## Default Windows Time Keeping
By default, Windows 11 uses the W32Time service to synchronize the system clock. For non-domain-joined computers, it polls `time.windows.com` every 1024 seconds (about 17 minutes) initially, adjusting up to every 32768 seconds (9 hours) based on clock stability. This provides loose accuracy suitable for general use like Kerberos authentication, but drifts can exceed seconds over time due to hardware clock inaccuracies and network variability, which is insufficient for precise astronomy timestamps. Domain-joined systems sync via the domain hierarchy instead.

## Configuring for Higher Accuracy
To achieve 100 ms or better accuracy (under normal wired network conditions), adjust W32Time settings via registry and commands as shown below. Focus on reducing poll intervals and using reliable NTP servers like `pool.ntp.org`. You will need administrator access to the computer to run these commands.

### Step 1: Set W32Time to Automatically Start Up
- Open Services (search "services.msc").  You will see tens or hundreds of services listed by "Name".
- Find the line with the name "Windows Time".  
- Double click on "Windows Time" to open the Windows Time Properties pop up settings.
- Set Startup type property to `Automatic`, and start the service if stopped.
- Click `Apply` and `OK` then close the pop up settings.
- Close the Services window by clicking on the `X` in the top right corner.

### Step 2: Configure NTP Servers and Sync Type

- Open the `Command Prompt` with the option `Run as Administrator`.
- Enter the following command: `w32tm /config /manualpeerlist:"tick.usnogps.navy.mil tock.usnogps.navy.mil pool.ntp.org" /syncfromflags:manual /reliable:yes /update`.  This sets multiple reliable NTP time servers (e.g., from ntp.org and the USNO).  I have found these time servers to be reliable but your experience may vary.  If your organization provides a stratum 1 time server, I would suggest using it in addition to the servers shown above.
- Restart service: `net stop w32time && net start w32time`
- Resync: `w32tm /resync`
- Close the `Command Prompt` window by clicking on the `X` in the top right corner.

### Step 3: Adjust the Windows Registry for Tighter Polling
- Open Registry Editor (regedit).
- Navigate to `HKLM\SYSTEM\CurrentControlSet\Services\W32Time\Config`
  - Set `MinPollInterval` to 6 (REG_DWORD, 2^6 = 64 seconds).
  - Set `MaxPollInterval` to 10 (REG_DWORD, 2^10 = 1024 seconds).
  - Set `UpdateInterval` to 57e40 (hexadecimal) or 360000 (decimal) (REG_DWORD).  This will cause the clock correction to be applied gradually over a one-hour time period.  The units for this setting are 10 ms increments.
  - Set `FrequencyCorrectRate` to 4 (REG_DWORD).
- Navigate to `HKLM\SYSTEM\CurrentControlSet\Services\W32Time\TimeProviders\NtpClient`
  - Set `SpecialPollInterval` to 400 (hexadecimal) or 1024 (decimal) (REG_DWORD).  The units for this setting are seconds.
- Apply: `w32tm /config /update`, then restart service and resync.
- These settings will cause the Windows time service to begin by polling the servers every 64 seconds.  This is as short as possible while avoiding having your computer placed on a black list by the servers for polling them too frequently.  As your computer's clock becomes more stable and well adjusted, the polling period will back off from 64 seconds to 128, 256, 512, and finally 1,024 seconds.  It will continue to evaluate its clock precision continuously during the observing session and adjust as necessary.  The Windows time service will choose one of the time servers to be it primary source of "truth" and will continue to monitor the entire set to detect any trouble with the chosen source.  The primary source will be chosen based on a small and highly repeatable network delay time.

### Step 4: Monitor and Verify
- Check status: `w32tm /query /status /verbose` (look for low phase offset).
- Test offset: `w32tm /stripchart /computer:pool.ntp.org /period:1024 /samples:10`.  The d value in the printed stripchart is always a positive value and indicates the estimated round trip network delay between your computer and the time server.  The o value may be positive or negative and indicates the estimated offset between your computer and the time server.  Positive values indicate your computer's clock is behind the time server's clock, while negative values indicate your computer's clock is ahead of the time server's clock.

## Tips from User Experiences
While I have never observed a meaningful offset between my computer and the time server since I've been using these settimgs, I automatically log the offset into my unified astro-log file using the `w32tm /query /status /verbose` command as the source of data for my log file just in case a question ever comes up.