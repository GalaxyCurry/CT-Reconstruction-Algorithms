CT-Reconstruction-Algorithms README
# CT-Reconstruction-Algorithms 🩻
This repository contains implementations of several classic CT image reconstruction algorithms. The goal is to provide clean, reproducible reference code for researchers and students to study, compare, and extend CT reconstruction methods.
## Directory Structure 📂
The repository is organized into three main folders, with detailed division as follows:
- **CBCT** 🛰️
  Contains various algorithms for cone-beam CT reconstruction. Currently, only the FDK algorithm is added, and other related algorithms are under planning and development.
- **Documents** 📄
  Stores academic papers related to all CT reconstruction algorithms included in the repository, providing theoretical support for algorithm implementation and research.
- **Filtered_Back_Projection(FBP)** 🔍
  Focuses on 2D image reconstruction, including 2D reconstruction algorithms for parallel beam and fan beam, corresponding to the related content in the Algorithms section.
## Algorithms 🧩
- Fan Beam 、Parallel beam 📡 
  - Filtered Back Projection (FBP) 🔍
- Cone Beam 🛰️
  - Feldkamp-Davis-Kress (FDK) 🧮
    - 3D cone-beam CT (CBCT) image reconstruction and core reconstruction pipeline
  - Iterative Reconstruction 🔄
    - 1
    - 2
    - 3
  - AI-based Reconstruction 🤖
    - 1
    - 2
## Code References 📚
This project is built upon and references the following open-source CT reconstruction toolboxes:
- **ASTRA Toolbox** 💎
  Advanced CPU/GPU accelerated library for projection and back-projection operations.
  GitHub: [astra-toolbox/astra-toolbox](https://github.com/astra-toolbox/astra-toolbox)
- **TIGRE Toolbox** ⭐
  Open-source cone-beam CT reconstruction library with FDK and iterative algorithms.
  GitHub: [CERN/TIGRE](https://github.com/CERN/TIGRE)
## Contributing 🤝
Contributions including new algorithms, bug fixes, and documentation improvements are welcome.
Please open an issue to discuss before submitting a pull request.
## Acknowledgments 🙏
Thanks to the original authors of the CT reconstruction algorithms.
This implementation is inspired and supported by open-source CT reconstruction frameworks.
## Contact 📩
If you have any questions, suggestions or cooperation intentions, please feel free to contact me via email:
lmy15933285944@163.com
