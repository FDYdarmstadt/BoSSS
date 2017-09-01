/* =======================================================================
Copyright 2017 Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.BoSSSpad {
    public static class DatabaseUtils {

        /// <summary>
        /// Create the BoSSS database directory structure in
        /// <see cref="targetDirectory"/> and registers the directory in the
        /// DBE config (<see cref="register_db.Program"/>). The structure is:
        ///  - `./data`
        ///  - `./timesteps`
        ///  - `./grids`
        ///  - `./sessions`
        /// </summary>
        /// <param name="dbDir">
        /// File system path to database; must be either non existent or an empty directory.
        /// </param>
        static public void CreateDatabase(string dbDir) {
            DirectoryInfo targetDirectory = new DirectoryInfo(dbDir);
            if (!targetDirectory.Exists) {
                targetDirectory.Create();
            } else {
                if (targetDirectory.GetFiles().Length > 0)
                    throw new ArgumentException("Must be empty.");
                if (targetDirectory.GetDirectories().Length > 0)
                    throw new ArgumentException("Must be empty.");
            }

            // Create structure
            Directory.CreateDirectory(Path.Combine(targetDirectory.FullName, "data"));
            Directory.CreateDirectory(Path.Combine(targetDirectory.FullName, "timesteps"));
            Directory.CreateDirectory(Path.Combine(targetDirectory.FullName, "grids"));
            Directory.CreateDirectory(Path.Combine(targetDirectory.FullName, "sessions"));
        }


        /// <summary>
        /// Checks if the given directory has the correct directory structure
        /// for a s BoSSS database.
        /// </summary>
        static public bool IsValidBoSSSDatabase(string dbDir) {
            DirectoryInfo databaseDirectory = new DirectoryInfo(dbDir);
            int score = 0;
            foreach (DirectoryInfo directory in databaseDirectory.GetDirectories()) {
                switch (directory.Name) {
                    case "data":
                    case "timesteps":
                    case "grids":
                    case "sessions":
                    score++;
                    break;
                }
            }

            return (score == 4);
        }

        /// <summary>
        /// Deletes a BoSSS database.
        /// </summary>
        static public void DeleteDatabase(string dbDir) {
            if(!IsValidBoSSSDatabase(dbDir)) {
                throw new ArgumentException("Not a valid BoSSS database.");
            }

            Directory.Delete(dbDir, true);
        }



    }
}
