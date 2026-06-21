
## 14.7 Finding Application files

If your application uses files that are not part of the saved program such as database files, configuration files, etc., the runtime version has to be able to locate these files. The [file_search_path/2](consulting.html#file_search_path/2) mechanism in combination with the **-p** `alias` command line argument provides a flexible mechanism for locating runtime files.
