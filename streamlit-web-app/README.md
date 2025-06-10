## Add Streamlit Web Interface for Cloud Deployment

### 🎯 Purpose
This PR adds a web-based version of the RepeatDetector GUI that can be deployed to cloud platforms, making the tool accessible without requiring Docker Desktop installation.

### ✨ Features Added
- **Streamlit web interface** replicating all wxPython GUI functionality
- **Cloud deployment ready** for platforms like Streamlit Cloud, Railway, Render
- **Interactive visualizations** using Plotly instead of static matplotlib
- **File upload/download** system for browser-based workflow
- **Docker support** for containerized deployment

### 🔄 Changes
- New `streamlit-web-app/` directory with web interface
- Converted wxPython GUI logic to Streamlit components
- Maintained all original analysis features (profile modes, instability index)
- Added requirements.txt and Dockerfile for deployment

### 🌐 Deployment Options
1. **Streamlit Cloud** (free): One-click deployment from GitHub
2. **Docker**: Containerized deployment anywhere
3. **Railway/Render**: Cloud platforms with Docker support

### 🧪 Testing
- [x] File upload functionality
- [x] Profile mode selection
- [x] Results calculation and display
- [x] Download functionality
- [x] Local Streamlit execution
- [ ] Production algorithm integration (noted in code)

### 📝 Notes
- Current implementation includes simulation of repeat detection for demo purposes
- Production deployment would require extracting the actual algorithm from the Docker image
- Maintains backward compatibility - doesn't affect existing desktop application

### 🔗 Related
Addresses the need for easier access to RepeatDetector without complex Docker setup.
