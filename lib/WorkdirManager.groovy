class WorkdirManager {

    static java.nio.file.Path resolveRunDir(def launchDir, String runName) {

        String today = java.time.LocalDate.now()
                        .format(java.time.format.DateTimeFormatter.ofPattern("yyyy-MM-dd"))

        java.nio.file.Path workBase = java.nio.file.Paths.get(launchDir.toString(), "work")

        // Count how many dated run folders already exist for today
        int runNumber = 1
        if (java.nio.file.Files.exists(workBase)) {
            long existing = java.nio.file.Files.list(workBase)
                .filter { java.nio.file.Files.isDirectory(it) }
                .filter { it.getFileName().toString().startsWith(today + "_run") }
                .count()
            runNumber = existing + 1
        }

        java.nio.file.Path candidate = workBase.resolve("${today}_run${runNumber}_${runName}")
        java.nio.file.Files.createDirectories(candidate)
        return candidate
    }
}