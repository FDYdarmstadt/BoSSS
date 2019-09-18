using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.FSI_Solver {
    public static class FSI_Tutorial {

        public static void Start() {
            Console.WriteLine("Hello! You started the FSI tutorial.");
            Console.WriteLine("Wir beginnen damit ein Kontrollobjekt für den Solver anzulegen.");
            Console.WriteLine("Darin werden alle Informationen gesammelt die für die Rechnung benötigt werden");
            Console.WriteLine("Dazu gehören die Definition des Rechengitters, die Eigenschaften von Fluid und Partikeln etc.");
            Console.WriteLine("Ein neues Kontrollobjekt wird wie folgt angelegt:");
            Console.WriteLine("FSI_Control(int degree, string projectName, string projectDescription, List<string> tags)");
            
        }
    }
}
