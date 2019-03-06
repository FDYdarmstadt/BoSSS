using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MPI.Wrappers;
using System.IO;
using Newtonsoft.Json;
using Newtonsoft.Json.Bson;

namespace BoSSS.Foundation.IO
{
    class VersionedSerializer : VectorDataSerializer
    {
        public VersionedSerializer(IFileSystemDriver driver) : base(driver)
        {

        }

    }
}
